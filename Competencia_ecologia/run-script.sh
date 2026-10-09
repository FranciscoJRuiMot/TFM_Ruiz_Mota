#!/bin/bash

# read pod identifier to point at the right input parameters:
pod_id=$((POD_ID + 1))

# create a folder for the output files in the persistent volume:
mkdir -p output/$pod_id

###############################################
### set parameters for building the network ###
###############################################
network_size=10
mean_connectance=0.5
stdev_connectance=0.2
mean_production=1
stdev_production=0.2
seed=$pod_id
network_file=$pod_id.csv

##############################################
### set parameters for avida configuration ###
##############################################
# set the species feeding on the external resource:
chemo_resources="A"
inflow=10
conversion_rate=0.75
# set the number of updates (i.e., time) until the end of the experiment:
exp_time=10000

#################################
### generate a random network ###
#################################
while true
do
    ./get_network.sh \
        "$network_size" "$mean_connectance" "$stdev_connectance" \
        "$mean_production" "$stdev_production" "$seed" "$network_file"

    # check if the network can be adapted to the format required for competition:
    adapted_to_competition=$(python check_adapted_competition.py "$network_file")

    # break only if output is non-empty and contains at least one valid row
    if [ -n "$adapted_to_competition" ] &&
       printf '%s\n' "$adapted_to_competition" |
       awk -v chemo_resources="$chemo_resources" '
           $1 != chemo_resources && $2 != chemo_resources { found=1; exit }
           END { exit !found }
       '
    then
        break
    fi

    # if not, (i.e., the network cannot be adapted to the competition format) increment seed:
    seed=$((seed + 1000))
done

##########################################################################
### prepare the network without competition for a competition scenario ###
##########################################################################
# get a species pair (excluding the one feeding on the external resource) to compete for the same resource:
pair=$(printf '%s\n' "$adapted_to_competition" | \
awk -v var="$chemo_resources" '$1 != var && $2 != var' | \
awk -v seed="$seed" 'BEGIN{srand(seed)} {print rand(), $0}' | \
sort -n | head -n 1 | cut -d' ' -f2-)

# create the resources_table.csv file using the adjacency matrix without competition:
read x y <<< "$pair"
awk -F',' -v x="$x" -v y="$y" '
BEGIN {
    OFS=","
    print "species_id,resource_id"
    id=1
}

NR > 1 {
    letter=$1

    if (letter==x || letter==y) {
        if (!(x in ids)) {
            ids[x]=id
            ids[y]=id
            id++
        }

        print letter, ids[letter]
    }
    else {
        print letter, id
        id++
    }
}
' $network_file > resources_table.csv

# convert the adjacency matrix without competition into the required one for competition:
awk -F',' -v pair="$pair" '
BEGIN{
    OFS=","
    split(pair, p, " ")
    keep = p[1]
    remove = p[2]
}

NR==1{
    # find column indices:
    for(i=1; i<=NF; i++){
        col[$i] = i
    }

    keep_col = col[keep]
    remove_col = col[remove]

    # print header without removed column:
    first=1
    for(i=1; i<=NF; i++){
        if(i != remove_col){
            if(!first) printf OFS
            printf "%s", $i
            first=0
        }
    }
    printf "\n"
    next
}

{
    # merge columns:
    $keep_col = sprintf("%.3f", $keep_col + $remove_col)

    # print row without removed column:
    first=1
    for(i=1; i<=NF; i++){
        if(i != remove_col){
            if(!first) printf OFS
            printf "%s", $i
            first=0
        }
    }
    printf "\n"
}
' $network_file > competition_matrix.csv

#####################################
### run avida without competition ###
#####################################
# run python script to configure and run avida:
python chemostat.py --adj_matrix=$network_file --chemo_resources=$chemo_resources --inflow=$inflow --conversion_rate=$conversion_rate --exp_time=$exp_time
# move output files to the persistent volume:
mv abundances.csv "output/${pod_id}/abundances_without_competition_${pod_id}.csv"
mv founders.csv "output/${pod_id}/founders_${pod_id}.csv"
mv "$network_file" "output/${pod_id}/network_without_competition_${pod_id}.csv"
# add also a new column indicating the pair of species that compete for the same resource:
awk -F',' -v pair="$pair" '
BEGIN{
    OFS=","
    split(pair, a, " ")
}
NR==1{
    print $0, "sp_competition"
    next
}
{
    flag=0
    for(i in a){
        if($1 == a[i]) flag=1
    }
    print $0, flag
}
' node_properties.csv > "output/${pod_id}/results_without_competition_${pod_id}.csv"
# remove results file without competition:
rm node_properties.csv
# remove data folder:
rm -R data

##################################
### run avida with competition ###
##################################
# set additional files:
mv competition_matrix.csv "${pod_id}.csv"
competition_matrix="${pod_id}.csv"
resources_table=resources_table.csv
founders="output/${pod_id}/founders_${pod_id}.csv"
# run python script:
python chemostat.py --adj_matrix=$competition_matrix --resources_table=$resources_table --organisms_path=$founders --chemo_resources=$chemo_resources --inflow=$inflow --conversion_rate=$conversion_rate --exp_time=$exp_time
# move output files to the persistent volume:
mv abundances.csv output/$pod_id/abundances_with_competition_$pod_id.csv
mv $competition_matrix output/$pod_id/network_with_competition_$pod_id.csv
mv resources_table.csv output/$pod_id/resources_table_competition_$pod_id.csv
# add also a new column indicating the pair of species that compete for the same resource:
awk -F',' -v pair="$pair" '
BEGIN{
    OFS=","
    split(pair, a, " ")
}
NR==1{
    print $0, "sp_competition"
    next
}
{
    flag=0
    for(i in a){
        if($1 == a[i]) flag=1
    }
    print $0, flag
}
' node_properties.csv > "output/${pod_id}/results_with_competition_${pod_id}.csv"
# remove results file with competition:
rm node_properties.csv
# remove data folder:
rm -R data
