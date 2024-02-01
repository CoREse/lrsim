#!/usr/bin/env python3
import sys

Version=0

Lengths=[]
IDs=[]#(length, ins count, del count)
if len(sys.argv)==1 or sys.argv[1]=="-":
    f=sys.stdin
else:
    f=open(sys.argv[1],"r")
Former=0
for line in f:
    if len(line)==0:
        continue
    sl=line.split()
    if len(sl)==3:
        if sl[0]=="RL":
            Length=int(sl[1])
            Count=int(sl[2])
            Lengths+=[(Length,Former)]*Count
            Former=Length
    elif len(sl)==4:
        if sl[0]=="ID":
            Length=int(sl[1])
            InsCount=int(sl[2])
            DelCount=int(sl[3])
            IDs+=[(Length,InsCount,DelCount)]
if len(sys.argv)==1 or sys.argv[1]=="-":
    f.close()
Ins=0
Del=0
for ID in IDs[:3]:
    Ins+=ID[0]*ID[1]
    Del+=ID[0]*ID[2]
Del/=Ins/100
Ins=100
Del=int(Del+0.5)
sys.stdout.buffer.write(Version.to_bytes(4,sys.byteorder))
sys.stdout.buffer.write(Ins.to_bytes(4,sys.byteorder))
sys.stdout.buffer.write(Del.to_bytes(4,sys.byteorder))
sys.stdout.buffer.write(len(Lengths).to_bytes(4,sys.byteorder))
for Length, Former in Lengths:
    sys.stdout.buffer.write(Length.to_bytes(4,sys.byteorder))
for Length, Former in Lengths:
    sys.stdout.buffer.write(Former.to_bytes(4,sys.byteorder))