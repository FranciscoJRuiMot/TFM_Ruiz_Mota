#!/bin/bash

# usage: ./get_network.sh N mean_C sd_C mean_prod sd_prod SEED output.csv
# N = number of species
# mean_C = mean connectance
# sd_C = standard deviation of connectance
# mean_prod = mean production per species
# sd_prod = standard deviation per species
# SEED = reproducibility
# output.csv = output file

if [ $# -ne 7 ]; then
    echo "usage: $0 <num_species> <mean_connectance> <sd_connectance> <mean_prod> <sd_prod> <seed> <output_file>"
    exit 1
fi

N=$1
MEAN_C=$2
SD_C=$3
MEAN_PROD=$4
SD_PROD=$5
SEED=$6
OUT=$7

# total possible directed links (excluding self-loops):
MAX_LINKS=$(( N*(N-1) ))

# generate connectance C from normal distribution with mean and SD:
C=$(awk -v mean=$MEAN_C -v sd=$SD_C -v seed=$SEED 'BEGIN{
    srand(seed)
    u1 = rand(); u2 = rand();
    z = sqrt(-2*log(u1)) * cos(2*3.14159*u2)
    c = mean + sd*z
    if (c < 0) c = 0
    if (c > 1) c = 1  # ensure C does not exceed 1
    printf "%.3f", c
}')

# calculate number of links based on C:
NUM_LINKS=$(awk -v max=$MAX_LINKS -v C=$C 'BEGIN{ printf("%d", max*C + 0.999) }')

# ensure minimum links to guarantee connectivity (random spanning tree):
CHAIN_LINKS=$((N-1))
if [ "$NUM_LINKS" -lt "$CHAIN_LINKS" ]; then
    NUM_LINKS=$CHAIN_LINKS
    C_actual=$(awk -v l=$CHAIN_LINKS -v max=$MAX_LINKS 'BEGIN{printf("%.3f", l/max)}')
    echo "sampled C too low; using minimum connectance C_actual = $C_actual"
else
    C_actual=$C
fi
#echo "$NUM_LINKS"
#echo "connectance = $C_actual"

# step 1: build a random spanning tree to guarantee connectivity:
> links.txt
for ((i=2;i<=N;i++)); do
    j=$((RANDOM % (i-1) + 1))
    echo "$j $i" >> links.txt
done
CHAIN_LINKS=$((N-1))

# step 2: remaining links randomly:
REMAINING_LINKS=$(( NUM_LINKS - CHAIN_LINKS ))
# generate all possible links (i, j), avoiding self-loops:
if [ $REMAINING_LINKS -gt 0 ]; then
    > possible_links.txt
    # create all possible links (i, j) where i != j:
    for ((i=1;i<=N;i++)); do
        for ((j=1;j<=N;j++)); do
            if [ $i -ne $j ]; then
                echo "$i $j" >> possible_links.txt
            fi
        done
    done
    # shuffle and select the remaining links (randomize):
    shuf --random-source=<(yes $SEED) possible_links.txt | while read link; do
        # add the link to links.txt if it has not been added already:
        if ! grep -q "$link" links.txt; then
            echo "$link" >> links.txt
            ((REMAINING_LINKS--))  # decrease remaining links to reach NUM_LINKS:
        fi
        # stop when we have reached the required number of links:
        if [ $REMAINING_LINKS -le 0 ]; then
            break
        fi
    done
fi

# step 3: create zero NxN matrix:
> tmp_matrix.txt
for ((i=1;i<=N;i++)); do
    row=$(yes 0 | head -n $N | tr '\n' ' ' | sed 's/ $//')
    echo $row >> tmp_matrix.txt
done

# step 4: assign random fractions per producer:
awk -v N=$N -v seed=$SEED -v mean_prod=$MEAN_PROD -v sd_prod=$SD_PROD '
BEGIN{
    # set the random seed to ensure reproducibility:
    srand(seed)
    
    # create species labels: A, B, C, ...
    for(i=1;i<=N;i++) 
        label[i]=sprintf("%c",64+i)
}
{
    # initialize matrix M to 0:
    for(i=1;i<=NF;i++) 
        M[NR,i]=0
}
FNR==NR{
    # read links (edges) and store them in "links" array:
    links[FNR]=$0
    next
}
END{
    # process links and build target species list for each species:
    for(l=1;l<=length(links);l++){
        split(links[l], arr)
        src=arr[1]
        tgt=arr[2]
        targets[src] = targets[src] " " tgt
    }

    # for each species, calculate production and proportions:
    for(i=1;i<=N;i++){
        split(targets[i], tarr)
        n=length(tarr)
        
        if(n==0) continue
        
        # generate species-specific production using Box-Muller method:
        u1 = rand()
        u2 = rand()
        z = sqrt(-2*log(u1)) * cos(2*3.14159*u2)
        total = mean_prod + sd_prod * z  # user-specified mean and std

        # avoid non-positive totals:
        if (total <= 0) 
            total = 0.01

        # generate random proportions for species (sum to 1):
        sumr = 0
        for(j=1;j<=n;j++) { 
            r[j] = rand()   # random number between 0 and 1
            sumr += r[j]    # sum of random numbers
        }

        # scale by total production:
        for(j=1;j<=n;j++){
            M[i,tarr[j]] = total * (r[j]/sumr)  # proportional scaling
        }
    }

    # print matrix with row & column labels
    # first row: empty cell + column labels
    header=""
    for(j=1;j<=N;j++){ header=header label[j] (j<N?",":"") }
    print "," header

    for(i=1;i<=N;i++){
        line=label[i] ","   # row label
        for(j=1;j<=N;j++){
            line=line sprintf("%.3f", M[i,j])
            if(j<N) line=line ","
        }
        print line
    }
}
' links.txt tmp_matrix.txt > $OUT

rm -f tmp_matrix.txt links.txt possible_links.txt
