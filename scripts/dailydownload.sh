#!/bin/bash

A="https://ncss.hycom.org/thredds/ncss/GLBv0.08/expt_53.X/data/2014?var=salinity&north=23&west=-159&east=-158&south=22&horizStride=1&time_start="
B="T00%3A00%3A00Z&time_end="
C="T23%3A59%3A59Z&timeStride=1&vertStride=1&accept=netcdf4"

# List of days in each month to handle the correct number of days
days_in_month=(31 28 31 30 31 30 31 31 30 31 30 31)

for month in {01..12}; do
  day_count=${days_in_month[$((10#$month - 1))]}
  for day in $(seq -w 01 $day_count); do
    date="2014-$month-$day"
    
    output="/Volumes/My Data/aloha/ALOHA Cabled Observatory Database/HYCOM/sal/2014/${month}_${day}.nc"
    
    # Construct the URL with the specific date for both time_start and time_end
    url="${A}${date}${B}${date}${C}"
    
    # Print the URL and output file for debugging
    echo "Downloading data for $date"
    echo "URL: $url"
    echo "Output file: $output"
    
    # Use curl to download the file
    curl "$url" -o "$output"
    
    # Check if the file was downloaded successfully
    if [ -f "$output" ]; then
      echo "Download completed for $date"
    else
      echo "Download failed for $date"
    fi
  done
done


##

#!/bin/bash

# Base URL components
A="https://ncss.hycom.org/thredds/ncss/GLBv0.08/expt_53.X/data/2014?var=water_temp&north=23&west=-159&east=-158&south=22&horizStride=1&time_start="
B="T00%3A00%3A00Z&time_end="
C="T23%3A59%3A59Z&timeStride=1&vertStride=1&accept=netcdf4"

# Directory to save the downloaded files
output_dir="/Volumes/My Data/aloha/ALOHA Cabled Observatory Database/HYCOM/tem/2014"

# Create the directory if it doesn't exist
mkdir -p "$output_dir"

# List of days in each month to handle the correct number of days
days_in_month=(31 28 31 30 31 30 31 31 30 31 30 31)

for month in {01..12}; do
  day_count=${days_in_month[$((10#$month - 1))]}
  for day in $(seq -w 01 $day_count); do
    date="2014-$month-$day"
    
    output="$output_dir/${month}_${day}.nc"
    
    # Construct the URL with the specific date for both time_start and time_end
    url="${A}${date}${B}${date}${C}"
    
    # Print the URL and output file for debugging
    echo "Downloading data for $date"
    echo "URL: $url"
    echo "Output file: $output"
    
    # Use curl to download the file
    curl "$url" -o "$output"
    
    # Check if the file was downloaded successfully
    if [ -f "$output" ]; then
      echo "Download completed for $date"
    else
      echo "Download failed for $date"
    fi
  done
done

##
#!/bin/bash

# Base URL components
A="https://ncss.hycom.org/thredds/ncss/GLBv0.08/expt_53.X/data/2014?var=surf_el&north=23&west=-159&east=-158&south=22&horizStride=1&time_start="
B="T00%3A00%3A00Z&time_end="
C="T23%3A59%3A59Z&timeStride=1&vertStride=1&accept=netcdf4"

# Directory to save the downloaded files
output_dir="/Volumes/My Data/aloha/ALOHA Cabled Observatory Database/HYCOM/ssh/2014"

# Create the directory if it doesn't exist
mkdir -p "$output_dir"

# List of days in each month to handle the correct number of days
days_in_month=(31 28 31 30 31 30 31 31 30 31 30 31)

for month in {01..12}; do
  day_count=${days_in_month[$((10#$month - 1))]}
  for day in $(seq -w 01 $day_count); do
    date="2014-$month-$day"
    
    output="$output_dir/${month}_${day}.nc"
    
    # Construct the URL with the specific date for both time_start and time_end
    url="${A}${date}${B}${date}${C}"
    
    # Print the URL and output file for debugging
    echo "Downloading data for $date"
    echo "URL: $url"
    echo "Output file: $output"
    
    # Use curl to download the file
    curl "$url" -o "$output"
    
    # Check if the file was downloaded successfully
    if [ -f "$output" ]; then
      echo "Download completed for $date"
    else
      echo "Download failed for $date"
    fi
  done
done
