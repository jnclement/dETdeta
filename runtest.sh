for i in range {0..200}; do
    root -b -q $1 2>&1 | grep -q "malloc"
    if [ $? == 0 ]; then
	echo "Malloc error found."
	break
    fi
    echo "Safe"
done
