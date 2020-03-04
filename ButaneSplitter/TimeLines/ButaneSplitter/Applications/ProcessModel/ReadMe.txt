Steps to create add on file with trend history on:
1. copy the model file (.mdl) to this folder
2. Verify that the correct model name is used in the Python-script MakeTrend.py
3. Update the K-SpiceTrendDef file to include all module-types variables you will add to trend with history
4. Run Python-script MakeTrend.py - this will read K-SpiceTrendDef.txt and Model-file, and create a new file: PP_AddOn.trn.
5. Copy the new file into your .trn file on your model-folder in between:

*********************

 </TrendDefinition>

"lines from new file in between here in the model.trn file"

 </TrendDefinitions>