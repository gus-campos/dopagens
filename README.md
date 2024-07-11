## O que é e para que serve?

Esse conjunto de programas serve para automatizar alguns processos envolvidos na otimização e avaliação de diversas moléculas estudadas pelo atual projeto de pesquisa.

Nesse projeto, usamos como base, 5 estruturas de grafino monocamada que têm de 1 a 5 ligações triplas entre os anéis benzênicos, e pretendemos substituir diversos átomos de carbono, de estrutura, por outros elementos, um a um, otimizar usando o DFTB+, e avalair os resultados.

<img src="IMAGENS/moleculas-base.jpg" alt="drawing" width="600"/>


Porém, para cada elemento, temos 70 versões difentes, o que tornaria o processo todo muito trabalhoso.

Portanto, esse programa tem o objetivo de automatizar alguns processos.

## Funções primárias

* Modificar dados xyz, em dado sítio, com dado elemento.

* Executar uma otimização, esperar ela acabar, e checar a convergência.

* De cada otimização, extrair informação sobre índex, elemento, coordenadas, e carga do átomo.

* Gera um histograma padronizado a partir do endereço de otimização.

* Retorna a energia final 

* Recebe endereço da otimização, chama função de energia total, e retorna energia de formação, de acordo com o elemento.

## Funções de lote

Usando as variáveis globais, que consideram a estrutura definida das dopagens, e considerando um elemento específico:

* Criar diretórios de dopagens.

* Modificar em lote os arquivos xyz (dopar), e escrevê-los no diretório de dopagem do elemento.

* Copiar arquivos de propriedades do dftb para os diretórios de dopagens do elemento.

* Iniciar uma fila de otimização com todas as moléculas de um elemento dopador.

* Gerar um arquivo csv da energias de formação de cada molécula dopada com tal elemento.

* Gerar arquivos padronizados de dados para todas as moléculas dopadas com esse elemento.

* A partir dos arquivo de dados, gerar histograma de ligações (e dat de histograma do grace) de cada molécula da dopagem com tal elemento.

* Usando os mesmos dados, gerar o mapa de carga de cada molécula das dopagens desse elemento.

## Estrutura de diretórios utilizada

```
.
├── arquivo
|
├── DOCS
│   └── IMAGENS
|   
├── DOP_INPUT
|   |
│   ├── CONFIG
|   |   |
│   │   ├── atoms_data.json
│   │   ├── config.json
│   │   └── glob_data.json
|   |
│   ├── COORDENADAS
│   ├── DFTB_IN
│   └── INDEXES
|   
├── ENERG_ATOM
|   |
│   └── [CLASSE]
|   
├── PROC_OUTPUT
|   |
│   ├── CHARGES_MAP
|   │   └── [CLASSE]
|   |   
│   ├── ENERGIAS
|   │   └── [CLASSE]
|   |   
│   ├── FRAMES
|   │   └── [CLASSE]
|   |   
│   ├── HISTOGRAMA
|   │   └── [CLASSE]
|   |   
│   └── HOMO_LUMO
|       └── [CLASSE]
|   
└── SCRIPTS

```

MAXANGULARMOMENTUM

B C N O F - “p”
Al Si P S Cl -> “d”
 