class ModelBlock:
        BlockType = ""
        TagName = ""
        def __init__(self):
                self.BlockType = ""
                self.TagName = ""
        def __init__(self,BlockTypee,TagNamer):
                self.BlockType = BlockTypee
                self.TagName = TagNamer
        def get_BlockType(self):
                return self.BlockType
        def set_BlockType(self, BlockTypee):
                self.BlockType = BlockTypee
        def get_TagName(self):
                return self.TagName
        def set_TagName(self, TagNamer):
                self.TagName = TagNamer
                
class HistoryVar:
        BlockType = ""
        BlockVar = ""
        def __init__(self):
                self.BlockType = ""
                self.BlockVar = ""
        def __init__(self,BlockTypee,BlockVarr):
                self.BlockType = BlockTypee
                self.BlockVar = BlockVarr
        def get_BlockType(self):
                return self.BlockType
        def set_BlockType(self, BlockTypee):
                self.BlockType = BlockTypee
        def get_BlockVar(self):
                return self.BlockVar
        def set_BlockVar(self, BlockVarr):
                self.BlockVar = BlockVarr

        
########################### 
#Array for holding all different types of history vars
HistoryVars = []

# Read File K-SpiceTrendDef.txt file
for line in open ("K-SpiceTrendDef.txt", "r").readlines():
        theLine = line.strip()
        theWords = theLine.split(":")
        HistoryVars.append(HistoryVar(theWords[0],theWords[1]))
        
###########################
# TEST

print(len(HistoryVars))
for h in HistoryVars:
#        print h.get_BlockType(), ":", h.get_BlockVar()
        

###########################
#Array for holding blocks in model
        ModelBlocks = []

# Read File K-Spice model file
for line in open ("ButaneSplitter.mdl", "r").readlines():
        l = line.strip()
        if l.find("<Block Name") > -1:
                block = l.split("\"")
                ModelBlocks.append(ModelBlock(block[3],block[1]))

###########################
# TEST          
#print len(ModelBlocks)

###########################
# - Scan the ModelBlocks
#       - Scan the HistoryVars
#               - if ModelBlocks[x].get_BlockType == HistoryVars[y].get_BlockType
#                       - add TrendItemDef for each HistoryVars[y] type
# - output to file
TrendDefs = ""
for m in ModelBlocks:
        for h in HistoryVars:
                if m.get_BlockType() == h.get_BlockType():
                        s = "<TrendItemDefinition name=\"" + m.get_TagName() + ":" + h.get_BlockVar() + "\" highScale=\"100\" lowScale=\"0\" zero=\"0\" description=\"\" />\n"
                        TrendDefs = TrendDefs + s

open ("PP_AddOn.trn", "w").write(TrendDefs)        

###########################
# End program
#print "TrendItemDefinitions now printed to TrendItemDefs.txt"
                        

# <TrendItemDefinition name="2FE1104:MeasuredValue" highScale="100" lowScale="0" zero="0" description="" />
        
  
