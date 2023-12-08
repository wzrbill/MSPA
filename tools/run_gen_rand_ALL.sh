echo "Out put k=5, |S|=5~8 and k=1e1~1e4, |S|=1e1~1e4."
for file in ./*; do
    if [[ "$file" == *".gr"* ]]; then
        ./rdgen $file 6     50
        echo "./rdgen $file 6     50"

        ./rdgen $file 7     50
        echo "./rdgen $file 7     50"

        ./rdgen $file 8     50
        echo "./rdgen $file 8     50"

        ./rdgen $file 9     50
        echo "./rdgen $file 9     50"

        ./rdgen $file 11     50
        echo "./rdgen $file 11     50"

        ./rdgen $file 101    50
        echo "./rdgen $file 101     50"

        ./rdgen $file 1001   50
        echo "./rdgen $file 1001     50"

        ./rdgen $file 10001  50
        echo "./rdgen $file 10001     50"
    fi
done

