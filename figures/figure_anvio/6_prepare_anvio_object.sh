#!/usr/bin/env bash
source activate anvio-env

# quit the session before adding more data:
anvi-import-misc-data viral_info.tsv -p anvio-object.db --target-data-table items --just-do-it
anvi-import-misc-data host_connection.tsv -p anvio-object.db --target-data-table items --just-do-it
anvi-import-misc-data host_info.tsv -p anvio-object.db --target-data-table layers --just-do-it
anvi-import-misc-data ordering_info.tsv -p anvio-object.db --target-data-table layer_orders --just-do-it
anvi-import-collection -p anvio-object.db -C vMAGs bins_collection.txt

