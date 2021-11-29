/bin/bash
mkdir -p "$HOME"/tmp/install
cd "$HOME"/tmp/install
wget https://github.com/brentp/hts-nim/releases/download/v0.2.8/hts_nim_static_builder && \
  chmod a+x hts_nim_static_builder

curl https://nim-lang.org/choosenim/init.sh -sSf > init.sh && \
  sh init.sh -y && \
  rm init.sh

export PATH=/root/.nimble/bin:$PATH

git clone https://github.com/brentp/hts-nim.git && \
  cd hts-nim && \
  nimble install -y
cd "$HOME"/tmp/install

git clone https://github.com/brentp/hts-nim-tools.git && \

hts_nim_static_builder -s hts-nim-tools/src/hts_nim_tools.nim -n hts-nim-tools/hts_nim_tools.nimble
