
outText = ""
tempWarn = "group = \"Deviation\"\n" + "set = 0\n" + "vtype = \"none\"\n" + "sample = 1\n" + "start = 0\n" + "alarm = 1\n" + "direct = 0\n" + "enum = \"\"\n"
tempAll = "ack = \"\"\n" + "enable = \"\";\n" + "actions = p1\n"

for line in open("alarmlist_dev.txt","r"):
        l = line.strip().split(";")
        tempH = "$dig \"" + l[0] + ":DeviationAlarm\"\n" + "synonym = \"" + l[1] + "\"\n" + "desc = \"High " + l[5] + "\"\n" + tempWarn + "set_ack = \"" + l[0] + ":ResetAlarmHigh\"\n" + tempAll
        outText = outText + tempH
        open ("AlarmDB_dev.dbl","w").write(outText)


print ("end")
 		
