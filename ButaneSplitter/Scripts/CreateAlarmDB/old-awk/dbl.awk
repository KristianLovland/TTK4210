BEGIN	{						# new file
		NR  = 1				# reset line nr.
		printf "###################################\n"
		printf "#                                 #\n"
		printf "#   ---------------------------   #\n"
		printf "#           Test.DBL              #\n"
		printf "#             FOR                 #\n"
		printf "#     Generic O&G PROJECT         #\n"
		printf "#   ---------------------------   #\n"
		printf "#                                 #\n"
		printf "#  ALARM TYPES:  1)  AlarmAnalog  #\n"
		printf "#                2)  SafetyValve  #\n"
		printf "#                                 #\n"
		printf "###################################\n"
		printf " \n"
	}
END   {
		printf " \n"
		printf "###DBL FILE FINISHED ###\n"
      }

NF > 0 {
		FS = "\t"				# make <TAB> the field separator
		gr0 = substr($3,1,2)  		# Read group addresse
		gr1 = substr($3,5,2)
		gr2 = substr($3,7,2)
		gr3 = substr($3,9,2)
#AlarmAnalog
		if ($4 == "AlarmAnalog")
		{
#	AlarmHighHigh
			printf("$dig %s:AlarmHighHigh\n", $1)
			printf("	synonym = \"%s\"\n",$2)
			printf("	desc = HH : \"%s\"\n",$6)
#			printf("	group = \"%s01%s%s%s\"\n", gr0,gr1,gr2,gr3)
			printf("	group = \"TripSignals\"\n")
			printf("   sample  = 1\n")
			printf("   alarm   = 1\n")
			printf("   direct  = 0\n")
			printf("   set_ack = \"%s:ResetAlarmHighHigh\"\n", $1)
			printf(";\n   actions = %s\n\n", $5)
#	AlarmHigh
			printf("$dig %s:AlarmHigh\n", $1)
			printf("	synonym = \"%s\"\n",$2)
			printf("	desc = H\n")
			printf("	group = \"Warnings\"\n")
			printf("   sample  = 1\n")
			printf("   alarm   = 1\n")
			printf("   direct  = 0\n")
			printf("   set_ack = \"%s:ResetAlarmHigh\"\n", $1)
			printf(";\n   actions = %s\n\n", $5)
#	AlarmLow
			printf("$dig %s:AlarmLow\n", $1)
			printf("	synonym = \"%s\"\n",$2)
			printf("	desc = L\n")
			printf("	group = \"Warnings\"\n")
			printf("   sample  = 1\n")
			printf("   alarm   = 1\n")
			printf("   direct  = 0\n")
			printf("   set_ack = \"%s:ResetAlarmLow\"\n", $1)
			printf(";\n   actions = %s\n\n", $5)
#	AlarmLowLow
			printf("$dig %s:AlarmLowLow\n", $1)
			printf("	synonym = \"%s\"\n",$2)
			printf("	desc = LL\n")
			printf("	group = \"TripSignals\"\n")
			printf("   sample  = 1\n")
			printf("   alarm   = 1\n")
			printf("   direct  = 0\n")
			printf("   set_ack = \"%s:ResetAlarmLowLow\"\n", $1)
			printf(";\n   actions = %s\n\n", $5)
		}
#SafetyValve
		if ($4 == "SafetyValve")
		{
#	open
			printf("$alm %s:ValveStemPosition\n", $1)
			printf("	synonym = \"%s\"\n",$1)
			printf("	desc = \"Open\"\n")
			printf("	group = \"%s01%s%s%s\"\n", gr0,gr1,gr2,gr3)
			printf("	set = 0\n")
			printf("	vtype = \"hi_hi\"\n")
			printf("	sample  = 1\n")
			printf("	enum = \"%s\"\n",$6)
			printf("	set_ack = \"%s\"\n",$6)
			printf("	ack = \"%s\"\n",$6)
			printf("	enable = \"%s\";\n",$6)
			printf("	hi_hi bell,p2,p1 = 0.01\n")
		}
	}
