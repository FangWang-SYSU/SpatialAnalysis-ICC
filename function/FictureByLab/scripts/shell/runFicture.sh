#!/bin/bash

#SBATCH --output=/home/%u/out/%x-%j.log
#SBATCH --time=80:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=7g
### Pixel level decoding
# Input: env, path, model_id, anchor_info, radius, coor, pixel

## param 1
mu_scale=1
major_axis="Y"
batch_size=500
batch_buff=30
skipshuffle=0


## param 2
key="Count"
train_width=10
sliding_step=2
min_ct_per_unit=20


## param 3
hex_width=12
anchor_res=4
nFactor=20
min_ct_per_feature=10
min_ct_per_unit_fit=10
cmap_name="turbo"
R=10
thread=8
train_nEpoch=10

# For DE output
n_top_gene=50
max_pval_output=1e-3
min_fold_output=1.5

pixel_resolution=0.8
pseudo='false'
pixel_resolution=0.6   # um, for visualization
analysis_resolution=0.1 # um

# For indexing pixel level output
bsize=2000
scale=100

plot_individual_factor=0
plot_subbatch=8


for ARGUMENT in "$@"
do
   KEY=$(echo $ARGUMENT | cut -f1 -d=)
   KEY_LENGTH=${#KEY}
   VALUE="${ARGUMENT:$KEY_LENGTH+1}"
   export "$KEY"="$VALUE"
done


set -xe
set -o pipefail

echo "Step.1 create batch file"
transcript=${in_path}/transcripts.tsv.gz
batchedOut=${in_path}/batched.matrix.tsv.gz
# pixel minibatch 
batch=$(echo $batchedOut | sed 's/\.gz$//g')
ficture make_spatial_minibatch --input ${transcript} --output ${batch} --mu_scale ${mu_scale} --batch_size ${batch_size} --batch_buff ${batch_buff} --major_axis ${major_axis}

if [ "${skipshuffle}" == "0" ]; then
   # shuffle minibatches
    sort -S 4G -k2,2g ${batch} | gzip -c > ${batch}.gz
    rm ${batch}
else
   gzip -f ${batch}
fi


# Prepare training minibatches, only need to run once if you plan to fit multiple models (say with different number of factors)
echo "Step 2. ficture make_dge: Create hexagon"
hexagon=${in_path}/hexagon.d_${train_width}.tsv.gz
out=$(echo $hexagon | sed 's/\.gz$//g')
ficture make_dge --key ${key} --count_header ${key} --input ${transcript} --output ${out} --hex_width ${train_width} --n_move ${sliding_step} --min_ct_per_unit ${min_ct_per_unit} --mu_scale ${mu_scale} --precision 2 --major_axis ${major_axis}

sort -S 4G -k1,1n ${out} | gzip -c > ${hexagon} # Shuffle hexagons
rm ${out}


echo "Step 3. Model training"
fit_nmove=${fit_nmove:-$((hex_width/anchor_res))}
model_id=${model_id:-nF${nFactor}.d_${train_width}}
output_id=${output_id:-${model_id}}
hexagon=${hexagon:-${in_path}/hexagon.d_${train_width}.tsv.gz}

anchor_info=prj_${hex_width}.r_${anchor_res}
output_path=${in_path}/analysis/${model_id}
figure_path=${output_path}/figure

if [ ! -d "${figure_path}/sub" ]; then
    mkdir -p ${figure_path}/sub
fi

output=${output_path}/${output_id}
#model=${output}.model.p

feature=${in_path}/feature.clean.tsv.gz

if [ "${pseudo}" = "true" ]; then
   # Fit model from pseudo-bulk data
   ficture init_model_from_pseudobulk --input ${hexagon} --output ${output}  --feature ${feature}  --model ${pseudo_matrix}  --epoch 0 --scale_model_rel -1  --key ${key}  --min_ct_per_feature ${min_ct_per_feature}  --thread ${thread} 
   cp ${output}.posterior.count.tsv.gz ${output}.model_matrix.tsv.gz
else
   # Fit model
   ficture fit_model --input ${hexagon}  --output ${output} --nFactor ${nFactor}  --feature ${feature}  --key ${key}  --epoch ${train_nEpoch}  \
                     --epoch_id_length 2  --test_split 0.5   --unit_attr X Y --min_ct_per_unit ${min_ct_per_unit}  --min_ct_per_feature ${min_ct_per_feature} --R ${R} --thread ${thread}
fi

# Choose color
input=${output_path}/${output_id}.fit_result.tsv.gz
output=${figure_path}/${output_id}


if [ -z "${cmap}" ]; then
    if [ -z "$color_anno" ]; then
        ficture choose_color --input ${input} --output ${output} --cmap_name ${cmap_name}
    else
        ficture choose_color --input ${input} --output ${output} --cmap_name ${cmap_name} --annotation ${color_anno}
    fi
else
    python /exdisk/zhangyanru/FICTURE/Script/plot_colors.py --input ${cmap} --output ${output}.png
fi


# Coarse plot for inspection
cmap=${cmap:-${figure_path}/${output_id}.rgb.tsv}
input=${output_path}/${output_id}.fit_result.tsv.gz
output=${figure_path}/${output_id}.coarse
fillr=$((hex_width/2+1))
ficture plot_base --input ${input} --output ${output} --fill_range ${fillr} --color_table ${cmap} --plot_um_per_pixel 1 --plot_discretized

# Transform
model=${model:-${output_path}/${output_id}.model_matrix.tsv.gz}
pixel=${transcript}
#pixel=${pixel:-${in_path}/filtered.matrix.tsv.gz}
output=${output_path}/${output_id}.${anchor_info}
ficture transform --input ${pixel} --output_pref ${output} --model ${model} --key ${key} --major_axis ${major_axis} --hex_width ${hex_width} --n_move ${fit_nmove} --min_ct_per_unit ${min_ct_per_unit_fit} --mu_scale ${mu_scale} --thread ${thread} --precision 2


 
echo "Step.4 "
coor=${coor:-${in_path}/coordinate_minmax.tsv}
pixel=${in_path}/batched.matrix.tsv.gz

output_path=${in_path}/analysis/${model_id}
figure_path=${output_path}/figure
output_id=${model_id}



anchor=${output_path}/${output_id}.${anchor_info}.fit_result.tsv.gz
K=$( echo $model_id | sed 's/nF\([0-9]\{1,\}\)\..*/\1/' )

#while IFS=$'\t' read -r r_key r_val; do
#    export "${r_key}"="${r_val}"
#done < ${coor}
#echo -e "${xmin}, ${xmax}; ${ymin}, ${ymax}"

while IFS=$'\t' read -r r_key r_val; do
    case $r_key in
        xmin) xmin=$(echo "${r_val} / ${mu_scale}" | bc -l) ;;
        xmax) xmax=$(echo "${r_val} / ${mu_scale}" | bc -l) ;;
        ymin) ymin=$(echo "${r_val} / ${mu_scale}" | bc -l) ;;
        ymax) ymax=$(echo "${r_val} / ${mu_scale}" | bc -l) ;;
        *) export "${r_key}"="${r_val}" ;;  # 其他键值对直接导出为环境变量
    esac
done < ${coor}
echo -e "${xmin}, ${xmax}; ${ymin}, ${ymax}"


# Projection without adaptation
radius=$(($anchor_res+1))
prefix=${output_id}.decode.${anchor_info}_${radius}
output=${output_path}/${prefix}
# Output only a few top factors per pixel
topk=3 # Fix for now
ficture slda_decode --input ${pixel} --output ${output} --model ${model} --feature ${feature} --anchor ${anchor} --anchor_in_um --neighbor_radius ${radius} --mu_scale ${mu_scale} --key ${key} --precision ${analysis_resolution}  --lite_topk_output_pixel ${topk} \
                    --lite_topk_output_anchor ${topk} --thread ${thread}


# Sort and index pixel level result
input=${output_path}/${prefix}.pixel.tsv.gz # j, X, Y, K1, ..., KJ, P1, ..., PJ, J=topk
output=${output_path}/${prefix}.pixel.sorted.tsv.gz

offsetx=${xmin}
offsety=${ymin}
rangex=$( echo "(${xmax} - ${xmin} + 0.5)/1+1" | bc )
rangey=$( echo "(${ymax} - ${ymin} + 0.5)/1+1" | bc )

header="##K=${K};TOPK=3\n##BLOCK_SIZE=${bsize};BLOCK_AXIS=X;INDEX_AXIS=Y\n##OFFSET_X=${offsetx};OFFSET_Y=${offsety};SIZE_X=${rangex};SIZE_Y=${rangey};SCALE=${scale}\n#BLOCK\tX\tY\tK1\tK2\tK3\tP1\tP2\tP3"

(echo -e "${header}" && zcat ${input} | tail -n +2 | perl -slane '$F[0]=int(($F[1]-$offx)/$bsize) * $bsize; $F[1]=int(($F[1]-$offx)*$scale); $F[1]=($F[1]>=0)?$F[1]:0; $F[2]=int(($F[2]-$offy)*$scale); $F[2]=($F[2]>=0)?$F[2]:0; print join("\t", @F);' -- -bsize=${bsize} -scale=${scale} -offx=${offsetx} -offy=${offsety} | sort -S 4G -k1,1g -k3,3g ) | bgzip -c > ${output}

tabix -f -s1 -b3 -e3 ${output}
rm ${input}

# DE
input=${output_path}/${prefix}.posterior.count.tsv.gz
output=${output_path}/DE/${prefix}.bulk_chisq.tsv
ficture de_bulk --input ${input} --output ${output} --min_ct_per_feature ${min_ct_per_feature} --max_pval_output ${max_pval_output} --min_fold_output ${min_fold_output} --thread ${thread}


# Report (color table and top DE genes)
#cmap=${output_path}/figure/${output_id}.rgb.tsv
#output=${output_path}/${prefix}.factor.info.html
if [ -z "$color_anno" ]; then
    ficture factor_report --path ${output_path} --pref ${prefix} --color_table ${cmap} --n_top_gene ${n_top_gene}
else
    ficture factor_report --path ${output_path} --pref ${prefix} --color_table ${cmap} --n_top_gene ${n_top_gene} --annotation ${color_anno}
    input=${output_path}/${prefix}.factor.info.tsv
    if [ ! -z ${ref_gene} ]; then
        python /exdisk/zhangyanru/FICTURE/Script/factor_report_change.py  --ficture_result ${input} --ref_gene ${ref_gene} --outpath ${output_path} --pref ${prefix} --column_name ${column_name}
    fi
fi


# Make pixel level figures
input=${output_path}/${prefix}.pixel.sorted.tsv.gz
output=${figure_path}/${prefix}.pixel.png
ficture plot_pixel_full --input ${input} --color_table ${cmap} --output ${output} --plot_um_per_pixel ${pixel_resolution} --full


if [ "${plot_individual_factor}" -eq 0 ]; then
    exit 0
fi

# If the data is very large, making all individual factor maps may take some time
# Make single factor heatmaps, plot_subbatch balances speed and memory ...
# batch size 8 should be safe for 7 or 14G in most cases
output=${figure_path}/sub/${prefix}.pixel.png
#if [ "${plot_subbatch}" -lt 1 ] || [ "${plot_subbatch}" -gt ${K} ]; then
#    ficture plot_pixel_single --input ${input} --output ${output} --plot_um_per_pixel ${pixel_resolution} --full --all
#else
#    st=0
#    ed=$((plot_subbatch+st-1))
#    while [ ${st} -lt ${K} ]; do
#        if [ ${ed} -gt ${K} ]; then
#            ed=$((K-1))
#        fi
#        id_list=$( seq ${st} ${ed} )
#        echo $id_list
#
#        ficture plot_pixel_single --input ${input} --output ${output} --id_list ${id_list} --plot_um_per_pixel ${pixel_resolution} --full
#
#        st=$((ed+1))
#        ed=$((plot_subbatch+st-1))
#    done
#fi

K=$(wc -l < ${cmap})
#id_list=$(seq 0 $((K-1)))
while IFS=$'\t' read -r name color_index r g b annotation; do
    if [[ "$name" == "Name" ]]; then
        continue
    fi
    r_int=$(echo "$r * 255" | bc)
    g_int=$(echo "$g * 255" | bc)
    b_int=$(echo "$b * 255" | bc)
    r_int=$(printf "%.0f" "$r_int")
    g_int=$(printf "%.0f" "$g_int")
    b_int=$(printf "%.0f" "$b_int")
    fixed_color="$r_int,$g_int,$b_int"
    echo "Fixed color for $name ($annotation): $fixed_color"
    current_id_list=$name
    ficture plot_pixel_single --input ${input} --output ${output} --id_list ${current_id_list} --plot_um_per_pixel ${pixel_resolution} --full --fixed_color ${fixed_color}
done < ${cmap}
