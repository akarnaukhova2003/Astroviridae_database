# Astroviridae_database
В папке data собранны данные по астровирусам. 
* interpro_input - в данной папке собраны 4 фаста-файла, которые подавались на вход interpro:
  - conflict_annot_ids.fasta (конфликтующие id)
  - not_annot_targetCDS_ids.fasta (id, для которых не нашлось таргетных рамок)
  - not_annot_CDS_nucl.fasta (нуклеотидные последовательности неаннотированных рамок)
  - not_annot_CDS_nucl_ORFs_after_getorfs.fasta (ак-последовательности,предсказанные с помощью get_orfs рамки для последовательностей из файла not_annot_CDS_nucl.fasta)
* interpro_выдача - в данной папке собраны файлы с выдачей interpro, а также результаты анализа данных файлов
  - archive - в этой папке собраны файлы, на основе которых мы уже записали координаты ORF, так что из них всю необходимую информацию получили. Не удаляю, чтобы при необходимости найти информацию по id, если вдруг будут проблемы с координатами
  - actual_for_191125 - папка с актуальными файлами
    - output_interpro_new - файла с сырой выдачей interpro: not_annot_targetCDS_new.tsv (для 48-ми геномов, для которых не нашлось таргетных ORF) и not_annot_CDS_nucl_ORFs_merged_output.csv (для 22-ти геномов, у которых не проаннотированы рамки)
    - predicted_orfs_after_interpro_new - тут собраны файлы, полученные в результате обработки выдачи интерпро: not_annot_CDS_nucl_ORFs_merged_output_orf_results_coords_AC*.csv и not_annot_targetCDS_orfs_new.csv - файлы, в которых я постаралась предсказать рамки из текстовой выдачи interpro; notannot_CDS_ORFs_findings.txt - собрала общую статистику по находнам для неаннотированных CDS
