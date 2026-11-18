#!/bin/bash

# 脚本开始时间
start_time=`date +%s`

# 创建有名管道并赋给文件描述符3
[ -e /tmp/fd1 ] || mkfifo /tmp/fd1
exec 3<>/tmp/fd1
rm -rf /tmp/fd1

# 输入用户参数
if [ ! -n "$1" ] ;then
    echo "You must input the chip id!";
    exit 1;
else
    chip=$1
fi

# 清空上次日志
#> $chip/report.log

#TODO: 改用getopt输入用户参数
for i in `seq 1 ${2:-32}`; do
    echo $i >&3 # 存入queue中，循环数限定为线程数量
done

# please move the raw data file under the folder './$id' into folder './$id/bcl'
bcl_in=$chip/bcl
fastq_in=$chip/fastq
mkdir $fastq_in
tilecoords_path=$chip/tiles_id
mkdir $tilecoords_path
file=$chip/bcl/RunInfo.xml
if [ ! -r "$file" ];
then
    chmod a+r $file
fi
if [ ! -w "$file" ];
then
    chmod a+w $file
fi

if [ -e $chip/lanes_and_tiles.txt ]; then
    echo $(date +"%Y-%m-%d %H:%M:%S") "[process: main]" 'Lanes and tiles coordinate file has existed!' >> $chip/report.log
else
    echo $(date +"%Y-%m-%d %H:%M:%S") "[process: main]" 'Creating lanes and tiles coordinate file!' >> $chip/report.log
    cat $file | grep "<Tile>" | sed 's/ *<Tile>//' | sed 's/<\/Tile>//' | sed 's/^[ \t]*//;s/[ \t]*$//' > $chip/lanes_and_tiles.txt
fi 

file=$chip/lanes_and_tiles.txt
if [ ! -r "$file" ];
then
    chmod a+r $file
fi
if [ ! -w "$file" ];
then
    chmod a+w $file
fi

# create tiles coordinate system

for tile in `cat $file`
do
    read -u3 thread
    suffix=fc_1_L${tile%_*}_tile_${tile##*_};
    echo $(date +"%Y-%m-%d %H:%M:%S") "[subprocess: $thread]" "Start: Demultiplex tile id $suffix" >> $chip/report.log
    {
        [ -d "$fastq_in/$suffix" ] || mkdir -p $fastq_in/$suffix;
        {
            bcl-convert --bcl-input-directory $bcl_in \
                     --output-directory $fastq_in/$suffix --tiles s_$tile \
                     --no-sample-sheet true -f
        } > $fastq_in/$suffix/report.log 2>&1
        echo $(date +"%Y-%m-%d %H:%M:%S") "[subprocess: $thread]" "Completed: Demultiplex tile id $suffix" >> $chip/report.log
        sleep 2
        echo $(date +"%Y-%m-%d %H:%M:%S") "[subprocess: $thread]" "Start: Map tile id $suffix coordinate system" >> $chip/report.log
        openst barcode_preprocessing \
            --fastq-in $fastq_in/$suffix/Undetermined_S0_L00${tile%_*}_R1_001.fastq.gz \
            --tilecoords-out $tilecoords_path \
            --out-suffix .txt \
            --out-prefix $suffix \
            --crop-seq 2:30 \
            --rev-comp \
            --single-tile >> $chip/report.log 2>&1
        sleep 2
        echo $(date +"%Y-%m-%d %H:%M:%S") "[subprocess: $thread]" "Completed: Map tile id $suffix coordinate system" >> $chip/report.log
        gzip $tilecoords_path/$suffix.txt >> $chip/report.log 2>&1
        echo $thread >&3
        echo $(date +"%Y-%m-%d %H:%M:%S") "[subprocess: $thread]" "subprocess-$thread: Completed!" >> $chip/report.log
    } &
    
done


# wait for subprocessing end
wait
stop_time=`date +%s`
time=$(( stop_time - start_time ))
echo $(date +"%Y-%m-%d %H:%M:%S") "[process: main]" "Runing time:${time}" >> $chip/report.log
exec 3<&-
exec 3>&-
echo flowcell tile barcode split finish! >> $chip/report.log