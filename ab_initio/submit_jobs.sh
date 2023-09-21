MAX_RUNNING=10000;
MAX_QUEUED=10;
SLEEP_TIME=1s;

function get_count_jobs {
    type=$1;
    echo `squeue -u $USER | grep " $type " | wc -l`;
}

function get_count_queued {
    get_count_jobs PD;
}

function get_count_running {
    get_count_jobs R;
}


while read jobcommand
do
    num_jobs_running=`get_count_running`;
    num_jobs_queued=`get_count_queued`;

    echo $num_jobs_running  $num_jobs_queued

    while (( $num_jobs_running >= $MAX_RUNNING  ||  $num_jobs_queued >= $MAX_QUEUED ));
    do
        sleep $SLEEP_TIME;
        num_jobs_running=`get_count_running`;
        num_jobs_queued=`get_count_queued`;
    done

    echo 'Submitting job';
    echo "> $jobcommand";
    eval $jobcommand;
done < "${1:-/dev/stdin}"
