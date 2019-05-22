# ProgramacaoDinamica
Código de Matlab que fornece políticas baseadas em Programação Dinâmica para o planejamento da operação energética de sistemas hidrotérmicos constituídos de um único reservatório. Também executa a simulação dessas políticas usando o histórico de vazões.

Dados de quatro usinas estão disponíveis no repositório: Furnas, Emborcação, Serra da Mesa e Sobradinho. Os dados podem ser alterados nos arquivos dados_furnas.m, dados_emborcacao.m, dados_serra_da_mesa.m e dados_sobradinho.m. 

O arquivo principal é o politicas.m, onde é possível escolher a usina para o estudo e o tipo de programação dinâmica (determinística, estocástica independente ou estocástica markoviana). Para mais detalhes sobre estas políticas, consulte a tese: https://www.researchgate.net/publication/333221626_Avaliacao_do_desempenho_de_politicas_em_regime_permanente_para_o_planejamento_da_operacao_energetica_de_sistemas_hidrotermicos_de_instancia_minima

Após sua execução, o algoritmo fornece os gráficos das políticas e das trajetórias de armazenamento obtidas por simulação. Também mostra as estatísticas da simulação sobre o histórico durante o período especificado.
