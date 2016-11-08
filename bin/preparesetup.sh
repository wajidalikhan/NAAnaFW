cmsenv
export XRD_NETWORKSTACK=IPv4
voms-proxy-init --voms cms
export FILECERT="x509up_u`id -u`"
echo $FILECERT
cp /tmp/$FILECERT $HOME
export X509_USER_PROXY=$HOME/$FILECERT
echo $X509_USER_PROXY
