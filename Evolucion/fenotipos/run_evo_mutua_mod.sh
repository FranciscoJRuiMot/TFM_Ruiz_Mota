#!/bin/bash

# Recibir el archivo CSV de entrada como argumento
input_csv=${1:-"abundances.csv"}
output_file="final_results"

# === CONFIGURACIÓN ===
env_id=1348
env_complexity="high_complexity"
input_1="0x0f0f5d44"
input_2="0x3366a3ad"
input_3="0x55728b55"
genome_size=100
# =====================

# Encabezado del archivo de resultados
echo "update,id,parent_id,abundance,genome_id,sequence,env_id,env_complexity,phenotype_id,NOT_timing,NOT_count,NAND_timing,NAND_count,AND_timing,AND_count,ORN_timing,ORN_count,OR_timing,OR_count,ANDN_timing,ANDN_count,NOR_timing,NOR_count,XOR_timing,XOR_count,EQU_timing,EQU_count,gest_time" > "$output_file.csv"

echo "=== Leyendo secuencias desde $input_csv ==="

# -----------------------------------------------------------------
# PASO 1: Leer TODAS las filas del CSV y construir:
#   - El fichero analyze.cfg con un LOAD_SEQUENCE por secuencia
#   - Arrays en memoria con los metadatos de cada fila
# -----------------------------------------------------------------
declare -a upd_arr
declare -a id_arr
declare -a parent_arr
declare -a abundance_arr
declare -a seq_arr

idx=0

# Construir la cabecera del .cfg
echo "PURGE_BATCH" > "analyze_${output_file}.cfg"

while IFS=',' read -r upd id parent abundance sequence; do
    echo "LOAD_SEQUENCE $sequence" >> "analyze_${output_file}.cfg"
    upd_arr[$idx]="$upd"
    id_arr[$idx]="$id"
    parent_arr[$idx]="$parent"
    abundance_arr[$idx]="$abundance"
    seq_arr[$idx]="$sequence"
    ((idx++))
done < <(tail -n +2 "$input_csv")

total=$idx
echo "=== $total secuencias cargadas en el .cfg ==="

# Añadir los comandos de análisis al .cfg (una sola vez, tras todos los LOAD_SEQUENCE)
echo "RECALC use_manual_inputs $input_1 $input_2 $input_3"                                                                                             >> "analyze_${output_file}.cfg"
echo "TRACE archive 0 -1 0 $input_1 $input_2 $input_3"                                                                                                 >> "analyze_${output_file}.cfg"
echo "DETAIL output.txt viable sequence gest_time task.8:binary task.7:binary task.6:binary task.5:binary task.4:binary task.3:binary task.2:binary task.1:binary task.0:binary" >> "analyze_${output_file}.cfg"

# -----------------------------------------------------------------
# PASO 2: Ejecutar Avida UNA SOLA VEZ con todas las secuencias
# -----------------------------------------------------------------
echo "=== Ejecutando Avida... ==="
./avida -a \
    -set DATA_DIR          "data_${output_file}" \
    -set ANALYZE_FILE      "analyze_${output_file}.cfg" \
    -set ENVIRONMENT_FILE  "environment_analyze.cfg" \
    -set VERBOSITY         0 \
    -set COPY_MUT_PROB     0 \
    -set DIVIDE_INS_PROB   0 \
    -set DIVIDE_DEL_PROB   0 \
    -set DIV_MUT_PROB      0 \
    -set OFFSPRING_SIZE_RANGE 1 \
    -set MIN_COPIED_LINES  0 \
    -set MIN_EXE_LINES     0 \
    -set REQUIRE_EXACT_COPY 1 \
    -set STERILIZE_UNSTABLE 1 \
    -set BASE_MERIT_METHOD 0 \
    -set BASE_CONST_MERIT  30 \
    -set DEATH_METHOD      0 \
    -set DIVIDE_RESET      1 > /dev/null 2>&1

echo "=== Avida terminado. Procesando trazas... ==="

# -----------------------------------------------------------------
# PASO 3: Iterar sobre cada secuencia usando su índice Avida
#         Avida nombra las trazas: org-Seq1.trace, org-Seq2.trace, ...
#         (índice base 1, por eso usamos seq_avida = idx+1)
# -----------------------------------------------------------------
for (( i=0; i<total; i++ )); do

    seq_avida=$(( i + 1 ))   # Avida numera desde 1
    trace_file="data_${output_file}/archive/org-Seq${seq_avida}.trace"

    upd="${upd_arr[$i]}"
    id="${id_arr[$i]}"
    parent="${parent_arr[$i]}"
    abundance="${abundance_arr[$i]}"
    sequence="${seq_arr[$i]}"

    # --- Viabilidad (línea 18 de output.txt, columna 1 == 1) ---
    # Nota: con múltiples secuencias output.txt tiene una línea por organismo;
    # la línea correspondiente a este organismo es seq_avida + 17 (cabecera de 17 líneas)
    # Si tu Avida escribe una línea por org en orden de carga, ajusta la lógica aquí.
    viable=$(awk -v line=$(( seq_avida + 17 )) 'NR==line && $1==1 {print $1}' "data_${output_file}/output.txt" 2>/dev/null)

    if [[ "$viable" -eq 1 ]]; then

        phen_decimal=$(awk -v line=$(( seq_avida + 17 )) '
            NR==line && $1==1 {
                bin=$4$5$6$7$8$9$10$11$12
                dec=0
                for(j=1; j<=length(bin); j++) dec = dec*2 + substr(bin,j,1)
                print dec
            }' "data_${output_file}/output.txt")

        gest_time=$(awk -v line=$(( seq_avida + 17 )) 'NR==line {print $3}' "data_${output_file}/output.txt")

        # --- Epigenoma desde la traza de ESTA secuencia ---
        awk '$2 ~ /IP:| / {print $2, $3}' "$trace_file" \
            | awk 'BEGIN{FS=":"} {print $1, $2, $3}' \
            | awk -v var=$genome_size '{if($2<var) print $3}' \
            | sed 's/(//' | sed 's/)//' \
            > "data_${output_file}/tmp_${seq_avida}.txt"

        epigenome_letters=$(awk 'BEGIN{FS=" "} NR==FNR{a[$1]=$2} NR>FNR{$1=a[$1];print}' \
            OFS='\n' instructions.txt "data_${output_file}/tmp_${seq_avida}.txt" | rs -Tg0)
        epigenome_length=${#epigenome_letters}

        if [[ "$gest_time" -eq "$epigenome_length" ]]; then

            # --- Timings: primera aparición de cada tarea en la traza ---
            NOT=$(awk  '/IP/ || /Task Count/ {print NR, $0}' "$trace_file" | awk '$5  == 1 && !found {print prev; found=1}{prev=$0}' | awk '{print $2}')
            NAND=$(awk '/IP/ || /Task Count/ {print NR, $0}' "$trace_file" | awk '$7  == 1 && !found {print prev; found=1}{prev=$0}' | awk '{print $2}')
            AND=$(awk  '/IP/ || /Task Count/ {print NR, $0}' "$trace_file" | awk '$9  == 1 && !found {print prev; found=1}{prev=$0}' | awk '{print $2}')
            ORN=$(awk  '/IP/ || /Task Count/ {print NR, $0}' "$trace_file" | awk '$11 == 1 && !found {print prev; found=1}{prev=$0}' | awk '{print $2}')
            OR=$(awk   '/IP/ || /Task Count/ {print NR, $0}' "$trace_file" | awk '$13 == 1 && !found {print prev; found=1}{prev=$0}' | awk '{print $2}')
            ANDN=$(awk '/IP/ || /Task Count/ {print NR, $0}' "$trace_file" | awk '$15 == 1 && !found {print prev; found=1}{prev=$0}' | awk '{print $2}')
            NOR=$(awk  '/IP/ || /Task Count/ {print NR, $0}' "$trace_file" | awk '$17 == 1 && !found {print prev; found=1}{prev=$0}' | awk '{print $2}')
            XOR=$(awk  '/IP/ || /Task Count/ {print NR, $0}' "$trace_file" | awk '$19 == 1 && !found {print prev; found=1}{prev=$0}' | awk '{print $2}')
            EQU=$(awk  '/IP/ || /Task Count/ {print NR, $0}' "$trace_file" | awk '$21 == 1 && !found {print prev; found=1}{prev=$0}' | awk '{print $2}')

            # --- Conteos finales de cada tarea ---
            NOT_count=$(grep  "Task Count (Quality):" "$trace_file" | tail -n 2 | head -n 1 | awk '{print $4}')
            NAND_count=$(grep "Task Count (Quality):" "$trace_file" | tail -n 2 | head -n 1 | awk '{print $6}')
            AND_count=$(grep  "Task Count (Quality):" "$trace_file" | tail -n 2 | head -n 1 | awk '{print $8}')
            ORN_count=$(grep  "Task Count (Quality):" "$trace_file" | tail -n 2 | head -n 1 | awk '{print $10}')
            OR_count=$(grep   "Task Count (Quality):" "$trace_file" | tail -n 2 | head -n 1 | awk '{print $12}')
            ANDN_count=$(grep "Task Count (Quality):" "$trace_file" | tail -n 2 | head -n 1 | awk '{print $14}')
            NOR_count=$(grep  "Task Count (Quality):" "$trace_file" | tail -n 2 | head -n 1 | awk '{print $16}')
            XOR_count=$(grep  "Task Count (Quality):" "$trace_file" | tail -n 2 | head -n 1 | awk '{print $18}')
            EQU_count=$(grep  "Task Count (Quality):" "$trace_file" | tail -n 2 | head -n 1 | awk '{print $20}')

            # --- Intercalar timing y count en una sola cadena ---
            values_line=""
            for var_name in NOT NAND AND ORN OR ANDN NOR XOR EQU; do
                var_value="${!var_name}"
                if [ -n "$var_value" ]; then values_line+="$var_value,"; else values_line+="NA,"; fi
            done
            values_line="${values_line%,}"

            values_line_count=""
            for count in "$NOT_count" "$NAND_count" "$AND_count" "$ORN_count" "$OR_count" "$ANDN_count" "$NOR_count" "$XOR_count" "$EQU_count"; do
                if [ -n "$count" ] && [ "$count" -ne 0 ] 2>/dev/null; then
                    values_line_count+="$count,"
                else
                    values_line_count+="NA,"
                fi
            done
            values_line_count="${values_line_count%,}"

            IFS=',' read -ra arr1 <<< "$values_line"
            IFS=',' read -ra arr2 <<< "$values_line_count"
            interleaved=""
            for j in "${!arr1[@]}"; do interleaved+="${arr1[$j]},${arr2[$j]},"; done
            interleaved="${interleaved%,}"

            # --- Escribir línea al CSV intermedio ---
            echo "$upd,$id,$parent,$abundance,$id,$sequence,$env_id,$env_complexity,$phen_decimal,$interleaved,$gest_time" >> "$output_file.csv"
        fi
    fi

    # Limpiar temporal de epigenoma de esta secuencia
    rm -f "data_${output_file}/tmp_${seq_avida}.txt"

done

# Limpiar directorio de datos y cfg tras procesar todo
rm -rf "data_${output_file}"
rm -f  "analyze_${output_file}.cfg"

# -----------------------------------------------------------------
# PASO FINAL: Generar columna 'phenotype' con AWK (igual que antes)
# -----------------------------------------------------------------
awk -F',' 'BEGIN{OFS=","}
NR==1 {print $0, "phenotype"; next}
{
    split("NOT,NAND,AND,ORN,OR,ANDN,NOR,XOR,EQU", tasks, ",")
    timing_cols[1]=10; timing_cols[2]=12; timing_cols[3]=14; timing_cols[4]=16
    timing_cols[5]=18; timing_cols[6]=20; timing_cols[7]=22; timing_cols[8]=24; timing_cols[9]=26
    count_cols[1]=11;  count_cols[2]=13;  count_cols[3]=15;  count_cols[4]=17
    count_cols[5]=19;  count_cols[6]=21;  count_cols[7]=23;  count_cols[8]=25;  count_cols[9]=27

    n_active = 0
    for (i=1; i<=9; i++) {
        t = $timing_cols[i]
        c = $count_cols[i]
        if (t != "NA" || c != "NA") {
            n_active++
            active_tasks[n_active] = tasks[i]
            active_timing[n_active] = (t != "NA") ? t : 999999
        }
    }
    for (i=1; i<=n_active; i++) {
        for (j=i+1; j<=n_active; j++) {
            if (active_timing[j] < active_timing[i]) {
                tmp_t = active_timing[i]; active_timing[i] = active_timing[j]; active_timing[j] = tmp_t
                tmp_n = active_tasks[i];  active_tasks[i] = active_tasks[j];  active_tasks[j] = tmp_n
            }
        }
    }
    phenotype = ""
    for (i=1; i<=n_active; i++) {
        phenotype = (i==1) ? active_tasks[i] : phenotype "+" active_tasks[i]
    }
    if (phenotype == "") phenotype = "NA"
    delete active_tasks
    delete active_timing
    print $0, phenotype
}' "$output_file.csv" > population_phenotype.csv

echo "=== Pipeline terminado. Resultado en population_phenotype.csv ==="