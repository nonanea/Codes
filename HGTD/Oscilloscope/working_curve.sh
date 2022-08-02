#!/bin/bash
file_dir=$1

rm $file_dir*txt
rm $file_dir*csv
txt_file=`find $file_dir -name "*results.txt"`

i=0
j=-1

# CT_array=(0 0 0 0)
# echo $txt_file
for obj in $txt_file
do
    echo $obj
    Bias=${obj%V*}
    Bias=${Bias##*_}
    # Bias=${obj%_*}
    # Bias=${Bias##*pix_}
    #echo $Bias"V"

    cat $obj | while read line
    do
        i=$(($i+1))

        if [[ $line = "amplitude_"* ]] || [[ $line = "charge_"* ]] || [[ $line = "RMS_"* ]] || [[ $line = "landau_"* ]] || [[ $line = "jitter_"* ]] || [[ $line = "sigma_"* ]] || [[ $line = "risetime_"* ]] || [[ $line = "tof_cfd_"* ]] || [[ $line = "SNR_"* ]]
        then
            # echo $i","$line
            
            j=$(($i+1))
            type=$line
        fi

        if [ $j = $i ]
        then
            is_mpv=$(echo $line | grep "MPV")
            if [ 0 = $? ]
            then
                mpv=${line#*MPV:}
                mpv=${mpv%%,*}
                
                mpv_err=${line#*MPV_ERR:}
                mpv_err=${mpv_err%%,*}

                # echo $Bias,$mpv,0,$mpv_err>>$file_dir$type.csv
                echo $Bias,$mpv,0,$mpv_err>>$file_dir$type.txt

                if [ $type = RMS_1 ]
                then
                    CT1=$mpv
                    CT3=$mpv_err
                fi

                if [ $type = charge_1 ]
                then
                    CT0=$mpv
                    CT2=$mpv_err

                    echo $CT0,$CT1,$CT2,$CT3>>$file_dir'NT'.txt
                fi
            fi

            is_mean=$(echo $line | grep "Mean")
            if [ 0 = $? ]
            then
                mean=${line#*Mean:}
                mean=${mean%%,*}
                
                mean_err=${line#*Mean_ERR:}
                mean_err=${mean_err%%,*}

                # echo $Bias,$mean,0,$mean_err>>$file_dir$type.csv
                echo $Bias,$mean,0,$mean_err>>$file_dir$type.txt
            fi

            if [ $type = "tof_cfd_23" ]
            then
                sigma=${line#*Sigma:}
                sigma=${sigma%%,*}

                sigma_err=${line#*Sigma_ERR:}
                sigma_err=${sigma_err%%,*}
                 
                # echo $Bias,$sigma,0,$sigma_err>>$file_dir$type-sigma.csv
                echo $Bias,$sigma,0,$sigma_err>>$file_dir$type-sigma.txt
            fi

            if [[ $type = "sigma_"* ]]
            then
                sigma=${line#*value:}
                sigma=${sigma%%,*}
                
                sigma_err=${line#*error:}
                sigma_err=${sigma_err%%,*}
                
                echo $Bias,$sigma,0,$sigma_err>>$file_dir$type.txt
                if [ $type = "sigma_2" ]
                    then
                        CT1=$sigma
                        CT3=$sigma_err

                        echo $CT0,$CT1,$CT2,$CT3>>$file_dir'CT'.txt
                        echo $Bias,$CT0,$CT1,$CT2,$CT3>>$file_dir'VCT'.txt
                fi
            fi

            # if [ $type = "sigma_1" -o $type = "sigma_2" -o $type = "sigma_3" ]
            if [[ $type = "jitter_ave_"* ]] || [[ $type = "landau_"* ]] || [[ $type = "risetime_ave"* ]] || [[ $type = "SNR_"* ]]
            then
                sigma=$line

                # echo $Bias,$sigma,0,0>>$file_dir$type.csv
                echo $Bias,$sigma,0,0>>$file_dir$type.txt
                    
                    if [ $type = SNR_1 ]
                    then
                        echo $CT0,$sigma,$CT2,0>>$file_dir'SNRC'.txt
                    fi
                    
                    if [ $type = jitter_ave_1 ]
                    then
                        echo $CT0,$sigma,$CT2,0>>$file_dir'JC'.txt
                    fi
            fi
        fi
    done
done