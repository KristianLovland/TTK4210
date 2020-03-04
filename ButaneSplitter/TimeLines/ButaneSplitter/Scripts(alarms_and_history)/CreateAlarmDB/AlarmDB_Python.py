
outText = ""
tempTrip = "group = \"Trip\"\n" + "set = 0\n" + "vtype = \"none\"\n" + "sample = 1\n" + "start = 0\n" + "alarm = 1\n" + "direct = 0\n" + "enum = \"\"\n"
tempWarn = "group = \"Warning\"\n" + "set = 0\n" + "vtype = \"none\"\n" + "sample = 1\n" + "start = 0\n" + "alarm = 1\n" + "direct = 0\n" + "enum = \"\"\n"
tempAll = "ack = \"\"\n" + "enable = \"\";\n" + "actions = p1\n"

for line in open("alarmlist.txt","r"):
        l = line.strip().split(";")
        tempHH = "$dig " + l[0] + ":AlarmHighHigh\"\n" + "synonym = \"" + l[1] + "\"\n" + "desc = \"HighHigh " + l[5] + "\n" + tempTrip + "set_ack = " + l[0] + ":ResetAlarmHighHigh\"\n" + tempAll
        tempH = "$dig " + l[0] + ":AlarmHigh\"\n" + "synonym = \"" + l[1] + "\"\n" + "desc = \"High " + l[5] + "\n" + tempWarn + "set_ack = " + l[0] + ":ResetAlarmHigh\"\n" + tempAll
        tempL = "$dig " + l[0] + ":AlarmLow\"\n" + "synonym = \"" + l[1] + "\"\n" + "desc = \"Low " + l[5] + "\n" + tempWarn + "set_ack = " + l[0] + ":ResetAlarmLow\"\n" + tempAll
        tempLL = "$dig " + l[0] + ":AlarmLowLow\"\n" + "synonym = \"" + l[1] + "\"\n" + "desc = \"LowLow " + l[5] + "\n" + tempTrip + "set_ack = " + l[0] + ":ResetAlarmLowLow\"\n" + tempAll
        outText = outText + tempHH + tempH + tempL + tempLL
        open ("AlarmDB.dbl","w").write(outText)


print ("end")
 		
