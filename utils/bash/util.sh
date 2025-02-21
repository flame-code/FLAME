inquire_yesno () {
    local answer=""
    while  ! ([ "$answer" == "y" ] || [ "$answer" == "n" ])
    do
        answer=""
        read -p "$1 (y/n)" answer
        answer=`echo "$answer" | tr '[:upper:]' '[:lower:]'`
        if [ "$answer" == "yes" ]; then answer="Y";fi
        if [ "$answer" == "no" ]; then answer="N";fi
    done
    if [ "$answer" == "y" ]; then
        echo true
    fi
    if [ "$answer" == "n" ]; then
        echo false
    fi
}
