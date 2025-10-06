struct Gate {
    (var id:Int,val type:String,var qb1:Int, var qb2:Int=-1, val param:String="") {
    var arity:Int
    var level=0
    var predecessors=ArrayList<Gate>()
    var successors=ArrayList<Gate>()

    init {
        arity=if(type=="SWAP" || qb2>=0) 2 else 1
    }

    override fun toString():String {
        val args = if(arity==1) "q[${qb1}]" else "q[${qb1}], q[${qb2}]"
        return "${type}${param} ${args}"
    }

};



class Problem {
    val gates=ArrayList<Gate>()
    var num_phases=1
    var depth=0
    var num_loc=0
    var num_res=0
    var num_gates=0
    var num_swaps=0
    val all_swaps=ArrayList<Pair<Int,Int>>()
    val measurements=ArrayList<Measure>()
    var preamble=""
    var uniform_cost=false
    lateinit var conn_matrix:IntArray
    lateinit var cost_matrix:IntArray
    lateinit var used_qubits:ArrayList<Int>
    lateinit var swap_ids:IntArray
    
    constructor(fname_arch:String, fname_gates:String, flag_read_precedences:Boolean=false) {
        load_architecture(fname_arch)
        //load_gates(fname_gates)
        load_qasm(fname_gates,flag_read_precedences)
        compute_cost_matrix()
    }
    
    constructor(p:Problem, gg:ArrayList<Gate>) {
        num_gates=0
        for(g in gg) {
            g.predecessors.clear()
            g.successors.clear()
            g.id=num_gates
            num_gates++
            gates.add(g)
        }
        num_loc=p.num_loc
        num_res=p.num_res
        for(sw in p.all_swaps) 
            all_swaps.add(sw)
        num_swaps=p.num_swaps
        for(m in p.measurements)
            measurements.add(m)
        preamble=p.preamble
        conn_matrix=p.conn_matrix
        cost_matrix=p.cost_matrix
        used_qubits=p.used_qubits
        swap_ids=p.swap_ids
        compute_precedences()
    }

    /*fun load_architecture_old(fname_arch:String) {
        val s=File(fname_arch).readText()
        val jo=JSONObject(s)
        num_loc=(jo["id"] as JSONObject).keySet().size
        val conn=jo["cx"] as JSONObject
        conn_matrix=IntArray(num_loc*num_loc)        
        for(c in conn.keys()) {
            val t=c.substring(1,c.length-1)
            val t1=t.split(",")
            val j1=t1[0].toInt()
            val j2=t1[1].toInt()
            conn_matrix[j1*num_loc+j2]=1
            conn_matrix[j2*num_loc+j1]=1
            all_swaps.add(Pair(j1,j2))
        }
        num_swaps=all_swaps.size
        swap_ids=IntArray(num_swaps*num_swaps, { -1 })
        var ix=0
        for(sw in all_swaps) {
            swap_ids[sw.first*num_swaps+sw.second]=ix
            swap_ids[sw.second*num_swaps+sw.first]=ix            
            ix++
        }
    }*/

    fun load_architecture(fname_arch:String) {
        val f=File(fname_arch)
        if(!f.exists()) {
            throw Throwable("${fname_arch} does not exist")
        }
        val ls=f.readLines()
        num_loc=ls[0].trim().toInt()
        for(i in 1..(ls.size-1)) {
    	    val ar=ls[i].trim().split(" ")
            val j1=ar[0].toInt()
            val j2=ar[1].toInt()
	        all_swaps.add(Pair(j1,j2))
        }
        num_swaps=all_swaps.size
        conn_matrix=IntArray(num_loc*num_loc, { 0 })        
        swap_ids=IntArray(num_swaps*num_swaps, { -1 })
        var ix=0
	    for((j1,j2) in all_swaps) {
            if(conn_matrix[j1*num_loc+j2]==0) {
                conn_matrix[j1*num_loc+j2]=1
                conn_matrix[j2*num_loc+j1]=1
                swap_ids[j1*num_swaps+j2]=ix
                swap_ids[j2*num_swaps+j1]=ix            
                ix++
            }
        }
        println("Loaded architecture with $num_loc qubits and $ix connections")
    }

    fun load_qasm(fname:String, flag_read_precedences:Boolean=false) {
        val f=File(fname)
        if(!f.exists()) {
            throw Throwable("${fname} does not exist")
        }
        val ls=f.readLines()
        used_qubits=ArrayList<Int>()
        var i=1
        val gate_regex=Regex("""\s*(\w+)\s*(\(.*\))?\s*q\[(\d+)\](\s*,\s*q\[(\d+)\])?""")
        val measure_regex=Regex("""\s*measure\s*q\[(\d+)\]\s*->c\[(\d+)\]""")
        num_gates=0
        var gates_unary=0
        var gates_binary=0
        while(i<ls.size) {
            val a=ls[i].trim().split(" ")
            if(a[0]=="//precedences")
                break
            if(a[0]=="OPENQASM" || a[0]=="include" || a[0]=="gate" || a[0]=="qreg" || a[0]=="creg")
                preamble+=ls[i]+"\n"
            else {
                val m1=measure_regex.find(ls[i])
                if(m1!=null) {
                    val gr=m1.groups
                    val qb=gr[1]!!.value.toInt()
                    val cb=gr[2]!!.value.toInt()                        
                    measurements.add(Measure(qb,cb))
                }
                else {
                    val m=gate_regex.find(ls[i])
                    if(m!=null) {
                        val gr=m.groups
                        //val id=a[0].toInt()
                        val type=gr[1]!!.value.uppercase()
                        val param=if(gr[2]!=null) gr[2]!!.value else ""
                        val qb1=gr[3]!!.value.toInt()
                        val qb2=if(gr[5]!=null) gr[5]!!.value.toInt() else -1
                        if(qb2>=0)
                            gates_binary++
                        else
                            gates_unary++
                        val g=Gate(num_gates,type,qb1,qb2,param)
                        gates.add(g)
                        num_gates++
                        if(!(qb1 in used_qubits)) used_qubits.add(qb1)
                        if(qb2>=0 && !(qb2 in used_qubits)) used_qubits.add(qb2)
                    }
                    else {
                        // errore
                        throw Throwable("I don't understand the line ${ls[i]}")
                    }
                }
            }
            i++
        }
        num_res=used_qubits.maxOrNull()!!+1
        //println("used qubits $used_qubits ($num_res)")
        if(i==ls.size || flag_read_precedences==false) 
            compute_precedences()
        else
            read_precedences(ls,i+1)
        println("Read circuit ${fname} with ${gates_unary} unary and ${gates_binary} binary gates, using ${num_res} qubits")
    }

    fun read_precedences(ls:List<String>,i0:Int) {
        println("reading precedences from qasm file")
        var i=i0
        val prec_regex=Regex("""//\s*(\d*)\s*(\d*)""")
        while(i<ls.size) {
            val res=prec_regex.find(ls[i])
            if(res!=null) {
                val i1=res.groups[1]!!.value.toInt()
                val i2=res.groups[2]!!.value.toInt()
                val g1=gates.find { it.id==i1 }
                val g2=gates.find { it.id==i2 }
                if(g1!=null && g2!=null) {
                    g1.successors.add(g2)
                    g2.predecessors.add(g1)
                }
            }
            i++
        }
        println("${i-i0} precedences read")
    }  

    fun compute_precedences() {
        println("computing precedences")
        val last=Array<Gate?>(num_res,{null})
        val tm=IntArray(num_res, {0})
        for(g in gates) {
            val p1=last[g.qb1]
            if(p1!=null) {
                p1.successors.add(g)
                g.predecessors.add(p1)
            }
            last[g.qb1]=g
            if(g.arity==2) {
                val p2=last[g.qb2]
                if(p2!=null) {
                    p2.successors.add(g)
                    g.predecessors.add(p2)
                }
                last[g.qb2]=g
                val st=1+Math.max(tm[g.qb1],tm[g.qb2])
                tm[g.qb1]=st; tm[g.qb2]=st
                g.level=st-1
            }
            else {
                g.level=tm[g.qb1]
                tm[g.qb1]++
            }
        }
        depth=tm.maxOrNull() ?: 0
        println("size $num_gates depth $depth")
    }

    fun load_gates(fname_gates:String) {
        val ls=File(fname_gates).readLines()
        used_qubits=ArrayList<Int>()
        var i=1
        while(true) {
            val a=ls[i].split(" ")
            if(a[0]=="precedences")
                break
            val id=a[0].toInt()
            //val phase=a[1].toInt()
            val type=a[1]
            val qb1=a[2].toInt()
            val qb2=a[3].toInt()
            //val param=a[4].toDouble()
            val g=Gate(id,type,qb1,qb2,a[4])
            gates.add(g)
            i++
            num_gates++
            if(!(qb1 in used_qubits)) used_qubits.add(qb1)
            if(qb2>=0 && !(qb2 in used_qubits)) used_qubits.add(qb2)
        }
        i++
        println("used qubits $used_qubits")
        num_res=used_qubits.size
        while(i<ls.size) {
            val a=ls[i].split(" ")
            val i1=a[0].toInt()
            val i2=a[1].toInt()
            val g1=gates.find { it.id==i1 }
            val g2=gates.find { it.id==i2 }
            if(g1!=null && g2!=null) {
                g1.successors.add(g2)
                g2.predecessors.add(g1)
            }
            i++
        }
    }
    
    fun adjacent(l1:Int, l2:Int) = 
        conn_matrix[l1*num_loc+l2]==1
    
    fun binary_gate_cost(l1:Int, l2:Int)=
        cost_matrix[l1*num_loc+l2]

    fun compute_cost_matrix() {
        fun set_cost(x:Int, y:Int, c:Int) {
            cost_matrix[x*num_loc+y]=c
        }
        cost_matrix=IntArray(num_loc*num_loc, {Int.MAX_VALUE})
        for (i in 0..num_loc-1) 
            //cost_matrix[i*num_loc+i]=0
            set_cost(i,i,0)
        for ((x,y) in all_swaps) {
            set_cost(x,y,1)
            set_cost(y,x,1)
        }
        var changed=true
        while (changed) {
            changed=false
            for (s in 0..num_loc-1)
                for (t in s+1..num_loc-1)
                    for ((x,y) in all_swaps) {
                        val d1=binary_gate_cost(s,x)
                        val d2=binary_gate_cost(y,t) 
                        if(d1<Int.MAX_VALUE && d2<Int.MAX_VALUE) {
                            val d=d1+d2+1
                            if(d<binary_gate_cost(s,t)) {
                                set_cost(s,t,d)
                                set_cost(t,s,d)
                                changed=true
                            }
                        }
                    }
        }
    }

    fun compute_duration(g:Gate,lo1:Int,lo2:Int=-1)=
        if(uniform_cost) 
            1
        else
            when(g.type) {
                "X","CX","H","MIX",
                "U1","U2","U3","T","CX","TDG",
                "RX","RY","RZ"    -> 1
                "WN_CX"           -> 4
                "PS"              -> 3
                "SWAP"            -> 1
                "MIXY"            -> 10 
                else  -> throw Throwable("Gate ${g.type} not valid")
            }


    fun get_swap_id(x:Int, y:Int):Int {
        val res=swap_ids[x*num_swaps+y]
        if (res==-1) throw Throwable("Swap on ($x,$y) is not possible")
        return res
    }
                    
}


