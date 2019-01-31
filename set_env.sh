export PROJDIR=$(git rev-parse --show-toplevel)
export R_ENVIRON_USER=${PROJDIR}/.Renviron
export R_PROFILE_USER=${PROJDIR}/.Rprofile

sed -i -e 's|R_LIBS=.*|R_LIBS=\"'"${PROJDIR}"'/.R_libs\"|' ${PROJDIR}/.Renviron
sed -i -e 's|R_LIBS_USER=.*|R_LIBS_USER=\"'"${PROJDIR}"'/.R_libs\"|' ${PROJDIR}/.Renviron
sed -i -e 's|DATADIR=.*|DATADIR=\"'"${PROJDIR}"'/data\"|' ${PROJDIR}/.Renviron
sed -i -e 's|RESDIR=.*|RESDIR=\"'"${PROJDIR}"'/result\"|' ${PROJDIR}/.Renviron

if [[ -n $NCORE ]]; then
    sed -i -e 's|NCORE=.*|NCORE='"${NCORE}"'|' ${PROJDIR}/.Renviron
fi