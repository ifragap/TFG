// Clear the output window
clear()

// Number of test and root to save data to set up - NEEDS TO CHANGE THE 
// ROOT IF NOT ON THE MAIN PC OF CTB

// If the "num_test" is 0, then a "cell_medium" folder will be created 
num_test = 1;             // THIS NEEDS TO CHANGE BETWEEN TESTS MANUALLY
test = "test" + num_test +"/";
const moment = new Date();

// First the time zone offset for later date export
offset = moment.getTimezoneOffset();
gmt = Math.abs(offset) / 60;
sign = offset > 0 ? '-' : '+';
gmt = (gmt < 10 ? "0" : "") + gmt;
minutes = Math.abs(offset) % 60;
minutes = (minutes < 10 ? "0" : "") + minutes;
off = "GMT" + sign + gmt + ":" + minutes;

// Then obtaining day, month and year
day = moment.getDate();
mth = moment.getMonth() + 1;
year = moment.getFullYear();
day = (day < 10 ? "0" : "") + day;
mth = (mth < 10 ? "0" : "") + mth;

// Creating part of the path of experiments
date = "" + day + "_" + mth + "_" + year + "/";

// These need to change if on another PC
path_1 = "D:/Inaki" + "/";
path_2 = "1/impedancia" + "/";

if(num_test == 0){
    root = path_1 + date + path_2 + "solo_medio/";
} else{
    root = path_1 + date + path_2 + test;
}

// Total time to measure is like: T_measure = X min * Y s * Z ms
time_min = 3;
time_seg = 60;
time_mil = 1000;
measure_time = time_min*time_seg*time_mil;

// Measures in X steps till we reach "measure_time"
for(var i = 0, end = Date.now() + measure_time; end > Date.now(); i++){
    t_start = Date.now();
    Impedance.single();
    t_end = Date.now();
    if(!Impedance.wait()) break;
    Impedance.Export(root + i + "_" + (t_start - (end - measure_time)) + "_" + (t_end - (end - measure_time)) + ".csv", "Impedance Analyzer", false);
}

// To save the date time we first get the start time of 
epoch = end - measure_time;
moment = new Date(epoch);

// Then the time
hr = moment.getHours();
min = moment.getMinutes();
sec = moment.getSeconds();
hr = (hr < 10 ? "0" : "") + hr;
min = (min < 10 ? "0" : "") + min;
sec = (sec < 10 ? "0" : "") + sec;

// Formation of the date + time with the pattern desired
date = off + " - " + day + "/" + mth + "/" + year + " - " + hr + ":" + min + ":" + sec + " - " + moment.getMilliseconds() + " ms";

// Converting the EPOCH time right
if(sign == '-'){
    epoch = epoch - (gmt*3600000 + min*60000);
} else{
    epoch = epoch + (gmt*3600000 + min*60000);
}

// Export of the file with the t_start in epoch and date
FileWrite(root + "t_start_test.txt", epoch + "\n" + date);