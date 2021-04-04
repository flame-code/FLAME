$1/$2.py $2.input $2.out.out 
awk 'BEGIN{i=0}{i++;printf("line%5.5d: %s\n",i,$0)}' $2.out.out >$2.out.yaml
sed -i 's/\t/ /g' $2.out.yaml
sed -i 's/  */ /g' $2.out.yaml
sed -i 's/ $/\r/g' $2.out.yaml
sed -i 's/$/]\r/g' $2.out.yaml
sed -i 's/ /,/g' $2.out.yaml
sed -i 's/:,/: [/g' $2.out.yaml
