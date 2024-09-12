# Changelog

## [0.3.0](https://github.com/RIVM-bioinformatics/juno-variant-typing/compare/v0.2.1...v0.3.0) (2024-09-12)


### Features

* add detection of snpCLs ([2d075d7](https://github.com/RIVM-bioinformatics/juno-variant-typing/commit/2d075d7a6bc8897425f3d28fc8047ccbff8e65bb))


### Bug Fixes

* add counts to lineage scoring ([d28718d](https://github.com/RIVM-bioinformatics/juno-variant-typing/commit/d28718d23fbdc2eb619eccd15f5cb7f771b160f0))
* dont use contig number in consensus fasta header ([514d806](https://github.com/RIVM-bioinformatics/juno-variant-typing/commit/514d8067666598cc9fce071f80e3cedc9b51197e))
* fix run_pipeline.sh ([e6dbab7](https://github.com/RIVM-bioinformatics/juno-variant-typing/commit/e6dbab75181cda89e780eb0428833526683a0755))
* fix run_pipeline.sh ([523a157](https://github.com/RIVM-bioinformatics/juno-variant-typing/commit/523a1571d99f15edcdd688241d6a8565414c7508))
* generalize binding paths when using container ([5d5991a](https://github.com/RIVM-bioinformatics/juno-variant-typing/commit/5d5991adaf45cc7a19a9bf9010170fae928b801f))
* introduce vcf entries with ref call if needed ([56ff964](https://github.com/RIVM-bioinformatics/juno-variant-typing/commit/56ff964ce005f21190177a2e9b30a520a3af44e6))


### Dependencies

* update juno-lib ([fbd67e6](https://github.com/RIVM-bioinformatics/juno-variant-typing/commit/fbd67e63bf67a7ae8f010c29a3df97626b43ba8f))

## [0.2.1](https://github.com/RIVM-bioinformatics/juno-variant-typing/compare/v0.2.0...v0.2.1) (2024-07-25)


### Bug Fixes

* update cli ([0de09d0](https://github.com/RIVM-bioinformatics/juno-variant-typing/commit/0de09d00d738929b61781f4bb8e38e23cafe981b))
* update cli ([57454c3](https://github.com/RIVM-bioinformatics/juno-variant-typing/commit/57454c375f678d758964893944c85c418308e69d))


### Dependencies

* update juno library ([97bb7ad](https://github.com/RIVM-bioinformatics/juno-variant-typing/commit/97bb7adbd6726e5e209922fbb7b4070dc9cdf44c))

## [0.2.0](https://github.com/RIVM-bioinformatics/juno-variant-typing/compare/v0.1.1...v0.2.0) (2024-07-05)


### Features

* annotate large deletions based on bed file ([095a2c7](https://github.com/RIVM-bioinformatics/juno-variant-typing/commit/095a2c7f85d450d17c75231463bb61cf85de1815))


### Bug Fixes

* dont fail when no deletions are found ([83d143f](https://github.com/RIVM-bioinformatics/juno-variant-typing/commit/83d143fcd472bb79be964fbfafb750a14cdb5ea6))

## [0.1.1](https://github.com/RIVM-bioinformatics/juno-variant-typing/compare/v0.1.0...v0.1.1) (2024-03-27)


### Bug Fixes

* allow multiple ref data rows merge for single variant ([8c4167e](https://github.com/RIVM-bioinformatics/juno-variant-typing/commit/8c4167e364f209d2712469329077a80007bff857))
* fix consensus generation ([f967cd5](https://github.com/RIVM-bioinformatics/juno-variant-typing/commit/f967cd5d804c29400dc8fda913f911cd8117406f))
* force alt regardless of GT ([a23a59b](https://github.com/RIVM-bioinformatics/juno-variant-typing/commit/a23a59b41f9ef75f8ca51a53c5ab77b69aa38b01))
* merge using python ([ec93c14](https://github.com/RIVM-bioinformatics/juno-variant-typing/commit/ec93c140ccf23de818dba6a2202c97a67ce5f336))
* set rules to local ([4a5e439](https://github.com/RIVM-bioinformatics/juno-variant-typing/commit/4a5e439f7203342318a87f6e03184de3d11bbeb5))

## 0.1.0 (2024-03-18)


### Features

* add consensus generation ([037d125](https://github.com/RIVM-bioinformatics/juno-variant-typing/commit/037d12503b4f30dab8380323827eeb5d1533e3ed))
* add manual version checks ([1577530](https://github.com/RIVM-bioinformatics/juno-variant-typing/commit/15775308830f71f4cfd40979541669869dd9bff3))
* add mtb lineage type as prototype ([f9f1273](https://github.com/RIVM-bioinformatics/juno-variant-typing/commit/f9f12734c88fafe961f479b4c53a36ac1eb204ea))
* add parsing of snpeff to variant table ([4cb2cbf](https://github.com/RIVM-bioinformatics/juno-variant-typing/commit/4cb2cbfad34f9922ecc4de9adf8056421799f572))
* add presets ([8c5309d](https://github.com/RIVM-bioinformatics/juno-variant-typing/commit/8c5309da7cb8c3d3a8a88c0e68d07f242be35c50))
* add small test dataset based on K1F phage ([d4a3bdb](https://github.com/RIVM-bioinformatics/juno-variant-typing/commit/d4a3bdbd0e468c0fa3a8eb1b6ca392a177d72ad0))
* add snpeff annotation ([a2123fa](https://github.com/RIVM-bioinformatics/juno-variant-typing/commit/a2123fa2e960bd36d95666d6249d0f53d732929f))
* copy and prepare input files ([4d98315](https://github.com/RIVM-bioinformatics/juno-variant-typing/commit/4d9831570368b1ccf1a5f5fa8d346c232bf9124f))
* include unfiltered variants in tsv output ([5832b92](https://github.com/RIVM-bioinformatics/juno-variant-typing/commit/5832b924543482ebce00e3e3b0625a89c287eb71))
* produce tb sequence experiment json for database ([502fc57](https://github.com/RIVM-bioinformatics/juno-variant-typing/commit/502fc57f6121795116ef8f281b99c5417f7b83cf))


### Bug Fixes

* add data to coll positions bed ([9765610](https://github.com/RIVM-bioinformatics/juno-variant-typing/commit/97656103e6ad09afa90dbb9ee719fa49e7fc4e04))
* add helper scripts for variant annotation ([88605d4](https://github.com/RIVM-bioinformatics/juno-variant-typing/commit/88605d42a80b9493ee52d474ded2083cf1425e93))
* add rrs and rrl count to output ([794a7a9](https://github.com/RIVM-bioinformatics/juno-variant-typing/commit/794a7a9ef25969e88a058e158f4053b71f1655df))
* bind db dir in wrapper ([df1c0c8](https://github.com/RIVM-bioinformatics/juno-variant-typing/commit/df1c0c81f790fdbb3487c2587e009e03105e8568))
* correct path to db ([2ac4337](https://github.com/RIVM-bioinformatics/juno-variant-typing/commit/2ac433798cfb19d2835273684dd3977c0d8e4e7c))
* correct path to resistance list on biogrid ([aea4ead](https://github.com/RIVM-bioinformatics/juno-variant-typing/commit/aea4ead9e82567857bc69ea1277f69b971e0dbfd))
* count only filtered variants in defined regions ([c7aa626](https://github.com/RIVM-bioinformatics/juno-variant-typing/commit/c7aa6268f316fd510bea8255e0c154f24857dba8))
* first conversion of template to jvt ([3b5d02d](https://github.com/RIVM-bioinformatics/juno-variant-typing/commit/3b5d02d91185aa4c84d7abb7999925c607ed6e60))
* fix error when no typing needed ([152bbc7](https://github.com/RIVM-bioinformatics/juno-variant-typing/commit/152bbc7bcce3d00c110195fa02e9711f6bd547d3))
* fix filtering for resistance mutations ([7d61cff](https://github.com/RIVM-bioinformatics/juno-variant-typing/commit/7d61cff7f6be8fa63ea90b23c89e3d9a9769f82c))
* fix mypy typing error ([e8f5118](https://github.com/RIVM-bioinformatics/juno-variant-typing/commit/e8f511838836cc7de2a2c94f5087c690fdaf9fe7))
* make column names more informative ([62ca364](https://github.com/RIVM-bioinformatics/juno-variant-typing/commit/62ca3640b750343865847b730dcac936befe6ef4))
* python executable name for container ([85a8ad2](https://github.com/RIVM-bioinformatics/juno-variant-typing/commit/85a8ad204d9f1fa3681387bd98653b947a6fd4cb))
* snpeff executable ([294eff7](https://github.com/RIVM-bioinformatics/juno-variant-typing/commit/294eff7eb9fe96dda18cf6188467b36741d39090))
* snpeff executable name ([2b9af66](https://github.com/RIVM-bioinformatics/juno-variant-typing/commit/2b9af669d2b71b5c0ad22d20c624890de35a4b66))
* sync conda and sing versions ([e82a948](https://github.com/RIVM-bioinformatics/juno-variant-typing/commit/e82a9487efa770d4dc5192178cc85f6a9be9020c))
* treat ref and allele columns as non-VCF compliant ([0e6a956](https://github.com/RIVM-bioinformatics/juno-variant-typing/commit/0e6a956b45db013b1cf954dd7e6ca706bff5f7ba))
* typo in default arg ([200210b](https://github.com/RIVM-bioinformatics/juno-variant-typing/commit/200210b03b02f25f40f2e5ed1415795879060b58))
* update typing ([06482cd](https://github.com/RIVM-bioinformatics/juno-variant-typing/commit/06482cd61df329cdb105f7617addb61f1e105020))


### Dependencies

* add conda envs ([5d2bb78](https://github.com/RIVM-bioinformatics/juno-variant-typing/commit/5d2bb78e8113b5f33a7eb748db9e8e350403b475))
* add containers ([d9d4f5e](https://github.com/RIVM-bioinformatics/juno-variant-typing/commit/d9d4f5e53ef38ef3214bc0e820e48db19f0a3a41))
* remove anaconda ([ca745ec](https://github.com/RIVM-bioinformatics/juno-variant-typing/commit/ca745eccc49d2c67b02af02bff411cfef989a5fc))
* set correct container ([cef23c4](https://github.com/RIVM-bioinformatics/juno-variant-typing/commit/cef23c413871466bee26ea0b1306022dc82783ad))
* sync conda install with containers ([f9a8d6b](https://github.com/RIVM-bioinformatics/juno-variant-typing/commit/f9a8d6b795ff8d1dab7e0ed14a6e64e0bf93a3a0))
* update library version ([8c23189](https://github.com/RIVM-bioinformatics/juno-variant-typing/commit/8c2318985437c5ef7ffd2b2de4f170f43ca9b6a2))
* upgrade snakemake and juno-lib ([cd94fed](https://github.com/RIVM-bioinformatics/juno-variant-typing/commit/cd94fed336197172a6555998d8ad41c237dc6168))
* use smk v7.24 to prevent reinstall with pip ([fca70fb](https://github.com/RIVM-bioinformatics/juno-variant-typing/commit/fca70fb59ed7be760cd098ebaefae71463bb18f9))
