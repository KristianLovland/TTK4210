
outText = ""
tempFail = "group = \"Failure\"\n" + "set = 0\n" + "vtype = \"none\"\n" + "sample = 1\n" + "start = 0\n" + "alarm = 1\n" + "direct = 0\n" + "enum = \"\"\n"
tempAll = "ack = \"\"\n" + "enable = \"\";\n" + "actions = p1\n"

for line in open("Alarms-Failures.txt","r"):
        l = line.strip().split(";")
        if l[3] == "ValveFailure":
                tempSF = "$dig \"" + l[0] + ":SignalFailure\"\n" + "synonym = \"" + l[1] + "\"\n" + "desc = \"Valve Signal Failure " + l[5] + "\"\n" + tempFail + "set_ack = \"\"\n" + tempAll
                tempF = "$dig \"" + l[0] + ":PowerFailure\"\n" + "synonym = \"" + l[1] + "\"\n" + "desc = \"Valve Power Failure " + l[5] + "\"\n" + tempFail + "set_ack = \"\"\n" + tempAll
        if l[3] == "MachineFailure":
                tempSF = "$dig \"" + l[0] + ":SignalFailure\"\n" + "synonym = \"" + l[1] + "\"\n" + "desc = \"Signal Failure " + l[5] + "\"\n" + tempFail + "set_ack = \"\"\n" + tempAll
                tempF = "$dig \"" + l[0] + ":MachineFailure\"\n" + "synonym = \"" + l[1] + "\"\n" + "desc = \"Machine Failure " + l[5] + "\"\n" + tempFail + "set_ack = \"\"\n" + tempAll
        outText = outText + tempSF + tempF
        open ("AlarmDB_Failures.dbl","w").write(outText)


print ("end")
 		
