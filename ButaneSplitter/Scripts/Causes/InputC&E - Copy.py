
class ChangeBlock:
 	tagName = ""
 	destText = ""
 	srcText = ""
 	def __init__(self,t,d,s):
 		self.tagName = t
 		self.destText = d
 		self.srcText = s
 	def getTagName (self):
 		return self.tagName
 	def getDestText (self):
 		return self.destText
 	def getSrcText (self):
 		return self.srcText

changeBlocks = []

for line in open("ModuleList.txt","r"):
	l = line.strip().split(";")
	c = ChangeBlock(l[0],l[1],l[2])
	changeBlocks.append(c)


print(len(changeBlocks))
print (changeBlocks[0].getTagName())
print (changeBlocks[0].getDestText())	
print (changeBlocks[0].getSrcText())

outText = ""
tagName = ""
appendInput = 0

inputText01 = "<Input>"
inputText02 = ""
inputText03 = ""
inputText04 = "</Input>"

for line in open("Model.mdl","r"):
        if appendInput > 1:
                if line.strip() == "<Inputs>":
                        tempText = line
                        print (tagName)
                                               
                        appendInput = 0
                        outText = outText + tempText + inputText01 + "\n" + destText + "\n" + srcText + "\n" + inputText04 + "\n"
                else:
                        outText = outText + line
        else:
                if line.count("<Block Name=\"") > 0:
                        tempText = line		
                        for block in changeBlocks:
                                if line.count("\""+block.getTagName()+"\"") > 0:
                                        l = line.split("\"")
                                        appendInput = 10
                                        tagName = block.getTagName()
                                        destText = block.getDestText()
                                        srcText = block.getSrcText()
                        outText = outText + tempText
                else:
                        outText = outText + line
 		
 
open ("Model_new.mdl","w").write(outText)
 		
