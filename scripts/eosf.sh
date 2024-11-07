#assumes lpc group space
theDir=lpclonglived
count_files_recursive() {
    local find_path=$1
    local ls_path=$2
    local count=0

    output=$(eos root://cmseos.fnal.gov find -f "$find_path" 2>&1)
    
    if echo "$output" | grep -q "Result is truncated"; then
        for subdir in $(eos root://cmseos.fnal.gov ls "$ls_path"); do
            sub_find_path="$find_path/$subdir"
            sub_ls_path="$ls_path/$subdir"
            if eos root://cmseos.fnal.gov ls "$sub_ls_path" > /dev/null 2>&1; then  # Check if it's a directory
                count=$((count + $(count_files_recursive "$sub_find_path" "$sub_ls_path")))
            fi  
        done
    else
        count=$(echo "$output" | wc -l)
    fi  

    echo $count
}

for dir in $(eos root://cmseos.fnal.gov ls /store/group/${theDir}); do
    top_find_path="/eos/uscms/store/user/${theDir}/$dir"
    top_ls_path="/store/group/${theDir}/$dir"
    total_count=$(count_files_recursive "$top_find_path" "$top_ls_path")
    echo "$dir $total_count files"
done
