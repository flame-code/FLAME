#true if the string $1 is inside $2
#to be controlled for empty strings
#stringContain() { [ -z "${2##*$1*}" ] && [ -z "$1" -o -n "$2" ]; }
stringContain() { test "${2##*$1*}" != "$2"; }

#define the value of the variable: unchanged if the string
#is already there otherwise prepended
getval()
{
    eval value=\$$1
    if stringContain $2 $value
    then
	val=
    else
	if test -z "$value"
	then
	    val=$2
	else
	    val=\"$2:\$"{$1}"\"
	fi
    fi
}

#portable command to set the environment. set $1 to $2
unienv()
{
    getval $1 $2
    if test -n "$val"
    then
	if stringContain 'bash' $SHELL
	then
	    toeval="export $1=$val"
	else
	    toeval="setenv $1 $val"
	fi
	echo $toeval ";"
	#eval $toeval #this only works in bash shells
    fi
}

