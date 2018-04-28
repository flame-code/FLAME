#for i in {1,5} 
for (( i=1; i <= 20; i++ ))
do
    vim -c sort -c wq fort.$((i+2000))
    vim -c sort -c wq fort.$((i+2100))
    vim -c sort -c wq fort.$((i+3000))
    vim -c sort -c wq fort.$((i+3100))
    echo "---------- diff fort.$((i+2000)) fort.$((i+2100))--------"
    diff fort.$((i+2000)) fort.$((i+2100))
    echo "---------- diff fort.$((i+3000)) fort.$((i+3100))--------"
    diff fort.$((i+3000)) fort.$((i+3100))

done
