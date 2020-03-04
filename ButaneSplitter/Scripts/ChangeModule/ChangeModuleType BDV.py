
class ChangeBlock:
 	blockType = ""
 	tagName = ""
 	blockDesc = ""
 	def __init__(self,b,t,d):
 		self.blockType = b
 		self.tagName = t
 		self.blockDesc = d
 	def getBlockType (self):
 		return self.blockType
 	def getTagName (self):
 		return self.tagName
 	def getBlockDesc (self):
 		return self.blockDesc

changeBlocks = []

for line in open("ModuleListBDV.txt","r"):
	l = line.strip().split(";")
	c = ChangeBlock(l[0],l[1],l[2])
	changeBlocks.append(c)

print(len(changeBlocks))
print (changeBlocks[2].getTagName())
	

outText = ""
 
for line in open("Model.mdl","r"):
 	if line.count("<Block Name=\"") > 0:
 		tempText = line
 		for block in changeBlocks:
 			if line.count("\""+block.getTagName()+"\"") > 0:
 				l = line.split("\"")
 				tempText = l[0] + "\"" + l[1] + "\"" + l[2] + "\"" + block.getBlockType() + "\"" + l[4] + "\"" + block.getBlockDesc() + "\"" + l[6] + "\"" + l[7] + "\"" + l[8] + "\""  + l[9] + "\"" + l[10]
 		outText = outText + tempText
 	else:
 		outText = outText + line
 		
 
open ("Model_new.mdl","w").write(outText)
 		
