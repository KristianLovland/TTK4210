outText = ""
for line in open("ModuleList.txt","r"):
    if line.strip() <> "":
        tempText = line
        outText = outText + line
open ("ModuleList_new.txt","w").write(outText)
