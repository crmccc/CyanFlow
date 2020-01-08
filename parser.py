# -*- coding:UTF8 -*-
#!/usr/bin/python3
import sys
import re
import getopt
class aline(object):
    def __init__(self,f,t,r,i):
        self.f=f
        self.t=t
        self.r=r
        self.i=i
class node:
    def __init__(self,node_type,a,b):
        '''0:PQ 1:PV 2:balance'''
        self.node_type=node_type
        self.a=a
        self.b=b

def main(argv):
    inputfile="alldatafile.txt"
    try:
        opts, args = getopt.getopt(argv,"hi:",["ifile="])
    except getopt.GetoptError:
        print 'parser.py <inputfile> '
        sys.exit()
    for opt, arg in opts:
        if opt == '-h':
            print 'parser.py <inputfile> '
            sys.exit()
        elif opt in ("-i", "--ifile"):
            inputfile = arg
    number_pattern=re.compile(r'-?\d+[.]\d+|\d')
    float_pattern=re.compile(r'-?\d+[.]\d+')
    complex_pattern=re.compile(r'-?\d+[.]\d+[+-]j\d+[.]\d+')
    file=open(inputfile,"r")
    lines=[]
    nodes=[]
    line_number=0
    pq_number=0
    pv_number=0
    number=0
    percision=0
    total=0
    for line in file:
        if line.find(r"题号")!=-1:
#output:
            if int(line.split('：')[1])==0:
                continue
            output=open("%d.txt" % number,'w')
            print(str(number))
            output.write("%d %d %d %d %f\n"% (total,pq_number,pv_number,line_number,percision))
            j=0
            for i in lines:
                j+=1
                output.write("%d %d %d %f %f\n" % (j,i.f,i.t,i.r,i.i))
            i=0
            for j in nodes:
                if j.node_type==0:
                    i+=1
                    output.write('%d %s %f %f\n'%(i,'Q',j.a,j.b))
                if j.node_type==1:
                    i+=1
                    output.write('%d %s %f %f\n'%(i,'V',j.a,j.b))
                if j.node_type==2:
                    i+=1
                    output.write('%d %s %f %f\n'%(i,'L',j.a,j.b))
            output.close()

            lines=[]
            nodes=[]
            line_number=0
            pq_number=0
            pv_number=0
            number=0
            total=0
            percision=0


            number=int(line.split('：')[1])

        if line.find(r"节点数：")!=-1:
            res=number_pattern.findall(line)
            total=int(res[0])
            line_number=int(res[1])
            percision=float(res[2])
        if line.find(r"支路 ")!=-1:
            res=complex_pattern.findall(line)
            res=str(res[0])
            res=res.replace('j','')
            res=res+'j'
            res=complex(res)
            r=res.real
            i=res.imag
        if line.find(r"─□─")!=-1:
            res=number_pattern.findall(line)
            f=int(res[0])
            t=int(res[1])
            lines.append(aline(f,t,r,i))
        if line.find(r"节点，")!=-1:
            res=number_pattern.findall(line)
            node_no=int(res[0])
            node_type=0
            s=0
            if line.find(r"ＰＱ节点，")!=-1:
                node_type=0
                res=complex_pattern.findall(line)
                res=str(res[0])
                res=res.replace('j','')
                res=res+'j'
                res=complex(res)
                pq_number+=1
                nodes.append(node(0,float(res.real),float(res.imag)))
            if line.find(r"ＰＶ节点，")!=-1:
                node_type=1
                res=float_pattern.findall(line)
                pv_number+=1
                nodes.append(node(1,float(res[0]),float(res[1])))
            if line.find(r"平衡节点，")!=-1:
                node_type=2
                res=float_pattern.findall(line)
                nodes.append(node(2,float(res[0]),float(res[1])))

    output=open("%d.txt" % number,'w')
    print(str(number))
    output.write("%d %d %d %d %f\n"% (total,pq_number,pv_number,line_number,percision))
    j=0
    for i in lines:
        j+=1
        output.write("%d %d %d %f %f\n" % (j,i.f,i.t,i.r,i.i))
    i=0
    for j in nodes:
        if j.node_type==0:
            i+=1
            output.write('%d %s %f %f\n'%(i,'Q',j.a,j.b))
        if j.node_type==1:
            i+=1
            output.write('%d %s %f %f\n'%(i,'V',j.a,j.b))
        if j.node_type==2:
            i+=1
            output.write('%d %s %f %f\n'%(i,'L',j.a,j.b))
    output.close() 

if __name__ == "__main__":
    main(sys.argv[1:])