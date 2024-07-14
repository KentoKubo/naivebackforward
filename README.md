# naivebackforward
naive forward/backward algorithm for pairHMM

## how to use

This program requires c++ environment.
I used docker container as a developing environment.

In your host system...
```bash
git clone git@github.com:KentoKubo/naivebackforward.git
docker-compose build
docker-compose run naivebackforward bash
```

In the container...
```bash
cd /workspace
make install
bf example/testseq.fa
```

- `bf` command caluculates forward/backward algorithm and shows the forward/backward score only. (now, only forward algorithm was inplemented)
- `bf` command requires a fasta file including 2 DNA sequences.
