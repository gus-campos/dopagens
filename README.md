# Dopings: Automação de dopagem, otimização e de análise de estruturas moleculares simuladas pelo DFTB+

Criei esse conjunto de programas, escritos em Python, como parte do meu projeto de pesquisa junto do Laboratório de Simulação Computacional do Departamento de Física da UFJF, com objetivo de pesquisar novos possíveis materiais e suas propriedades teóricas ao simular computacionalmente flocos de grafeno e grafino com impurezas diversas, em variadas posições. As simulações foram feitas utilizando o pacote de simulação quântica DFTB+.

Posteriormente ainda foram desenvolvido outros programas que geram versões de tais estruturas, adsorvidas com moléculas de H2, para estudar a possibilidade de usar tais estruturas em tecnologia de armazenamento de hidrogênio.

| Elementos | Al,  B, Li, Mg,  N,  Na,  O,  P, Si,  Ti,  Zn |
|-|-|
| Bases grafino | grafino 1, 2, 3, 4, 5 |
|Bases grafeno  | armchair edges, zizag edges |

Foram geradas, otimizadas e processadas, um total de 1.540 estruturas. Além de mais de 3.500 versões de estruturas adsorvidas com H2.

## Funcionalidades

* Criar versões em lote, de várias estruturas base, dopadas por vários elementos químicos, em diversas posições;
* Criar um diretório de otimização para cada estrutura, com os arquivos e informações necessárias;
* Criar filas de otimizações do DFTB+, com um leve nível de auto gerência;
* Gerar relatório do estado e progresso das otimizações;
* Extrair os dados necessário dos arquivos;
* Gerar árvores de diretórios organizadas com o a geometria final da otimização;
* Gerar visualizações e processamentos dos dados para viabilizar a análise da grande quantidade de resultados;
* Gerar novas estruturas adsorvidas com H2, para diversos tamanhos e formas de estruturas, tanto para floco, quanto para periódica.

## Dopagens

| ![](assets/g1.png) | ![](assets/N-g1-B2.png) | 
| - | - |
| ![](assets/ac.png) | ![](assets/P-ac-A1.png) |

## Adsorção com H2

| ![](assets/g5-150.png) | ![](assets/ac-200-v.png) | 
| - | - |
| ![](assets/g1_s4.png) | ![](assets/g1_s4-100.png) |

## Relatório de otimização

```
../opt_files/graphine/Zn/g5/C4 - Converged      | Last step: 277

        Not converged SCC's: [0, 1, 2, 3, 4, 5, 6, 7, 8] => 3.2% of total

../opt_files/graphine/Zn/g5/D1 - Converged      | Last step: 489
../opt_files/graphine/Zn/g5/D2 - Converged      | Last step: 395
../opt_files/graphine/Zn/g5/D3 - Converged      | Last step: 492
../opt_files/graphine/Zn/g5/D4 - Converged      | Last step: 397

        Not converged SCC's: [2, 3, 4] => 0.8% of total

../opt_files/graphine/Zn/g5/D5 - Converged      | Last step: 424

        Not converged SCC's: [2] => 0.2% of total

../opt_files/graphene/Al/ac/A1 - Converged      | Last step: 320
../opt_files/graphene/Al/ac/A2 - Converged      | Last step: 321
../opt_files/graphene/Al/ac/A3 - Converged      | Last step: 319
../opt_files/graphene/Al/ac/B1 - Converged      | Last step: 314
../opt_files/graphene/Al/ac/B2 - Converged      | Last step: 308
../opt_files/graphene/Al/ac/B3 - Converged      | Last step: 300
../opt_files/graphene/Al/ac/B4 - Converged      | Last step: 302
../opt_files/graphene/Al/ac/B5 - Converged      | Last step: 295
../opt_files/graphene/Al/ac/B6 - Converged      | Last step: 309
../opt_files/graphene/Al/ac/B7 - Converged      | Last step: 305
../opt_files/graphene/Al/ac/B8 - Converged      | Last step: 306
```

## Resultados

### Visualizações

| ![](assets/histogram.png) | ![](assets/charges_map.png) | 
| - | - |
| ![](assets/gap_graph.png) | ![](assets/geometry_graph.png) |




