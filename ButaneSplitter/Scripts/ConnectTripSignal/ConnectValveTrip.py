
class ChangeBlock:
 	tagName = ""
 	def __init__(self,t):
 		self.tagName = t
 	def getTagName (self):
 		return self.tagName

changeBlocks = []

for line in open("ModuleList.txt","r"):
	l = line.strip().split(";")
	c = ChangeBlock(l[0])
	changeBlocks.append(c)

print(len(changeBlocks))
print (changeBlocks[0].getTagName())
	

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
                        inputText02 = "<Destination>" + tagName + ":LocalSetClosed</Destination> \n"
                        inputText03 = "<Source>" + tagName + ":Tripped</Source> \n"
                        appendInput = 0
                        outText = outText + tempText + inputText01 + "\n" + inputText02 + inputText03 + inputText04 + "\n"
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
                        outText = outText + tempText
                else:
                        outText = outText + line
 		
 
open ("Model_new.mdl","w").write(outText)
 		
