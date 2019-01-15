#!/bin/sh

GREP_OPTIONS=''

cookiejar=$(mktemp cookies.XXXXXXXXXX)
netrc=$(mktemp netrc.XXXXXXXXXX)
chmod 0600 "$cookiejar" "$netrc"
finish() {
  rm -rf "$cookiejar" "$netrc"
}

trap finish EXIT
WGETRC="$wgetrc"

prompt_credentials() {
    echo "Enter your Earthdata Login or other provider supplied credentials"
    read -p "Username (ptaconet): " username
    username=${username:-ptaconet}
    read -p "Password: " password
    echo "machine urs.earthdata.nasa.gov login $username password $password" >> $netrc
    echo
}

exit_with_error() {
    echo
    echo "Unable to Retrieve Data"
    echo
    echo $1
    echo
    echo "https://e4ftl01.cr.usgs.gov//MODV6_Cmp_A/MOLA/MYD11A2.006/2018.08.29/MYD11A2.A2018241.h17v07.006.2018250040311.hdf"
    echo
    exit 1
}

prompt_credentials
  detect_app_approval() {
    approved=`curl -s -b "$cookiejar" -c "$cookiejar" -L --max-redirs 2 --netrc-file "$netrc" https://e4ftl01.cr.usgs.gov//MODV6_Cmp_A/MOLA/MYD11A2.006/2018.08.29/MYD11A2.A2018241.h17v07.006.2018250040311.hdf -w %{http_code} | tail  -1`
    if [ "$approved" -ne "302" ]; then
        # User didn't approve the app. Direct users to approve the app in URS
        exit_with_error "Please ensure that you have authorized the remote application by visiting the link below "
    fi
}

setup_auth_curl() {
    # Firstly, check if it require URS authentication
    status=$(curl -s -z "$(date)" -w %{http_code} https://e4ftl01.cr.usgs.gov//MODV6_Cmp_A/MOLA/MYD11A2.006/2018.08.29/MYD11A2.A2018241.h17v07.006.2018250040311.hdf | tail -1)
    if [[ "$status" -ne "200" && "$status" -ne "304" ]]; then
        # URS authentication is required. Now further check if the application/remote service is approved.
        detect_app_approval
    fi
}

setup_auth_wget() {
    # The safest way to auth via curl is netrc. Note: there's no checking or feedback
    # if login is unsuccessful
    touch ~/.netrc
    chmod 0600 ~/.netrc
    credentials=$(grep 'machine urs.earthdata.nasa.gov' ~/.netrc)
    if [ -z "$credentials" ]; then
        cat "$netrc" >> ~/.netrc
    fi
}

    fetch_urls() {
    if command -v curl >/dev/null 2>&1; then
        setup_auth_curl
        while read -r line; do
            curl -f -b "$cookiejar" -c "$cookiejar" -L --netrc-file "$netrc" -Og -- $line && echo || exit_with_error "Command failed with error. Please retrieve the data manually."
        done;
    elif command -v wget >/dev/null 2>&1; then
        # We can't use wget to poke provider server to get info whether or not URS was integrated without download at least one of the files.
        echo
        echo "WARNING: Can't find curl, use wget instead."
        echo "WARNING: Script may not correctly identify Earthdata Login integrations."
        echo
        setup_auth_wget
        while read -r line; do
        wget --load-cookies "$cookiejar" --save-cookies "$cookiejar" --keep-session-cookies -- $line && echo || exit_with_error "Command failed with error. Please retrieve the data manually."
        done;
    else
        exit_with_error "Error: Could not find a command-line downloader.  Please install curl or wget"
    fi
}

fetch_urls <<'EDSCEOF'
http://e4ftl01.cr.usgs.gov//MODV6_Dal_D/SRTM/SRTMGL1.003/2000.02.11/N11W004.SRTMGL1.hgt.zip
http://e4ftl01.cr.usgs.gov//MODV6_Dal_D/SRTM/SRTMGL1.003/2000.02.11/N10W004.SRTMGL1.hgt.zip
http://e4ftl01.cr.usgs.gov//MODV6_Dal_D/SRTM/SRTMGL1.003/2000.02.11/N09W006.SRTMGL1.hgt.zip
http://e4ftl01.cr.usgs.gov//MODV6_Dal_D/SRTM/SRTMGL1.003/2000.02.11/N08W006.SRTMGL1.hgt.zip
EDSCEOF
