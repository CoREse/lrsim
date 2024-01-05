import sys

Lengths=[]
if len(sys.argv)==1 or sys.argv[1]=="-":
    f=sys.stdin
else:
    f=open(sys.argv[1],"r")
Former=0
for line in f:
    sl=line.split()
    Length=int(sl[0])
    Count=int(sl[1])
    Lengths+=[(Length,Former)]*Count
    Former=Length
if len(sys.argv)==1 or sys.argv[1]=="-":
    f.close()
sys.stdout.buffer.write(len(Lengths).to_bytes(4,sys.byteorder))
for Length, Former in Lengths:
    sys.stdout.buffer.write(Length.to_bytes(4,sys.byteorder))
for Length, Former in Lengths:
    sys.stdout.buffer.write(Former.to_bytes(4,sys.byteorder))