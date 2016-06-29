#!/bin/bash

############################################################################################
# 35.runMAPS_InChrArmChunks.sh		Modified from 35.runMAPS_Pipeline_parallel_run_pileup.sh		
#
# Paul Bailey	5.5.2015
#
############################################################################################

referenceFile=/tgac/workarea/group-cg/baileyp/IWGSC_v2_ChrU_Ref/IWGSC_v2_ChrU.fa					# The full and final reference for the project (IWGSC CSS minus CSS 3B + new 3BSeq + CSS 3B not in 3BSeq + ChrU)
																									# Not used in MAPS chunking mode
splitRefBedfileLocatn=/tgac/workarea/group-cg/baileyp/references/wheat_genomes/my_iwgsc_v2/my_iwgsc_v2_split	# Location of the split reference 																							
	
#referenceFile=/tgac/workarea/collaborators/dubcovskylab/reference/IWGSC_CSS_AB-TGAC_UCW_v1.fa		# Reference for Ksenia and Hans' Kronos samples
#splitRefBedfileLocatn=/tgac/workarea/group-cg/baileyp/references/wheat_genomes/IWGSC_CSS_AB-TGAC_UCW_v1_fromKsenia/IWGSC_CSS_AB-TGAC_UCW_v1_split	# Location of the split reference 


#chrArms=(1AL 1AS 1BL 1BS 1DL 1DS 2AL 2AS 2BL 2BS 2DL 2DS 3AL 3AS 3B 3DL 3DS 4AL 4AS 4BL 4BS 4DL 4DS 5AL 5AS 5BL 5BS 5DL 5DS 6AL 6AS 6BL 6BS 6DL 6DS 7AL 7AS 7BL 7BS 7DL 7DS TGAC_Cadenza_U) # 41 arms + chrU
chrArms=(5DS 6AL)
#chrArms=(1AL 1AS 1BL 1BS 2AL 2AS 2BL 2BS 3AL 3AS 3B 4AL 4AS 4BL 4BS 5AL 5AS 5BL 5BS 6AL 6AS 6BL 6BS 7AL 7AS 7BL 7BS UCW_Kronos_U) 	# Reference for Ksenia and Hans' Kronos samples


# Produce the bam and bai files for each set of chromosome arm contigs (RUN TIME: 2h; a 27 sample data set took 7h to run; max 4h for mm_sampe_sort.bam):
:<<***string_marker_for_commenting_out_code***

# Location of the files: 
#bamFileLoctn=/tgac/workarea/group-cg/baileyp/WheatLoLa/35.ExonCapture_myIWGSC_v2_chrU_ref/54_test_samples
#bamFileLoctn=/tgac/workarea/group-cg/baileyp/WheatLoLa/35.ExonCapture_myIWGSC_v2_chrU_ref/Plate[A]
#	NB - I changed the name of some PlateA Sample_1147 samples - so need to use *sampe_sort_markdup_rm.bam for these.
#bamFileLoctn=/tgac/scratch/baileyp/WheatLoLa/35.ExonCapture_myIWGSC_v2_chrU_ref/PlateD
#bamFileLoctn=/tgac/scratch/baileyp/WheatLoLa/39.ExonCapture_PlateH2N/PlateH
bamFileLoctn=/tgac/workarea/group-cg/baileyp/WheatLoLa/39.ExonCapture_PlateH2N/PlateK
#bamFileLoctn=/tgac/workarea/group-cg/baileyp/WheatLoLa/40.MergingKronosData4Jorge

bamFileName=sampe_sort_markdup_rm.bam
#bamFileName=mm_sampe_sort.bam


# List of files:
#fileList=/tgac/workarea/group-cg/baileyp/WheatLoLa/35.ExonCapture_myIWGSC_v2_chrU_ref/54_test_samples/35.ExonCapture_PreparingTable_54_TestSamples.txt
#fileList=/tgac/workarea/group-cg/baileyp/WheatLoLa/35.ExonCapture_300_Samples/35.ExonCapture_PreparingTable.txt
fileList=/tgac/workarea/group-cg/baileyp/WheatLoLa/39.ExonCapture_PlateH2N/ExonCapture_PlateH2N.txt
#fileList=/tgac/workarea/group-cg/baileyp/WheatLoLa/40.MergingKronosData4Jorge/kronos_file_list.txt

#plateToDo='1128_LIB9960_LDI8215\|1128_LIB9961_LDI8216\|1128_LIB9962_LDI8217\|1128_LIB9963_LDI8218\|1157_LIB10000_LDI8255\|1157_LIB10001_LDI8256\|1157_LIB10002_LDI8257\|1157_LIB10003_LDI8258'
#plateToDo='1147_LIB'
#plateToDo='.'
# 54 test samples:
#plateToDo='LIB5772\|LIB5773\|LIB5774\|LIB5775\|LIB5776\|LIB5777\|LIB5778\|LIB5779\|LIB5787\|LIB5788\|LIB5789\|LIB5790\|LIB8410\|LIB8411\|LIB8413\|LIB8415\|LIB8419\|LIB8420\|LIB8422\|LIB8423\|LIB8424\|LIB8425\|LIB8426\|LIB8427'
#plateToDo='LIB8428\|LIB8429\|LIB8430\|LIB8431\|LIB8432\|LIB8433\|LIB8434\|LIB8435\|LIB8438\|LIB8439\|LIB8440\|LIB9205\|LIB9206\|LIB9207\|LIB9209\|LIB9210\|LIB9211\|LIB9212\|LIB9221\|LIB9222\|LIB9225\|LIB9226\|LIB9227\|LIB9231\|LIB9232'
#plateToDo='LIB9208\|LIB9228'	#'LIB5780\|LIB5782\|LIB8418'
# Plate A:
#plateToDo='LIB10604\|LIB10605\|LIB10606\|LIB10608\|LIB10609\|LIB10610\|LIB10611\|LIB10616\|LIB10617\|LIB10618\|LIB10619\|LIB10620\|LIB10621\|LIB10622\|LIB10623\|LIB10624\|LIB10625\|LIB10626\|LIB10627\|LIB10632\|LIB10633\|LIB10634\|LIB10635\|LIB10636\|LIB10638\|LIB10701'
#plateToDo='LIB10630\|LIB10640\|LIB10641\|LIB10642\|LIB10643\|LIB10656\|LIB10658\|LIB10659\|LIB10660\|LIB10661\|LIB10662\|LIB10663\|LIB10680\|LIB10681\|LIB10682\|LIB10683\|LIB10688\|LIB10690\|LIB10691\|LIB10692\|LIB10693\|LIB10696\|LIB10698\|LIB10699\|LIB10702\|LIB10705\|LIB10706'
# Plate B:
#plateToDo='LIB10020\|LIB10021\|LIB9924\|LIB9925\|LIB9926\|LIB9927\|LIB9928\|LIB9930\|LIB9931\|LIB9932\|LIB9933\|LIB9935\|LIB9936\|LIB9938\|LIB9939\|LIB9948\|LIB9949\|LIB9950\|LIB9951\|LIB9952\|LIB9953\|LIB9954\|LIB9955\|LIB9960\|LIB9961\|LIB9962\|LIB9963'
#plateToDo='LIB9964\|LIB9965\|LIB9966\|LIB9967\|LIB9968\|LIB9969\|LIB9970\|LIB9971\|LIB9972\|LIB9973\|LIB9974\|LIB9975\|LIB9976\|LIB9977\|LIB9978\|LIB9979\|LIB9984\|LIB9985\|LIB9986\|LIB9987\|LIB9988\|LIB9989\|LIB9990\|LIB9991\|LIB9992\|LIB9993\|LIB9994\|LIB9995'
#plateToDo='LIB10000\|LIB10001\|LIB10002\|LIB10003\|LIB10004\|LIB10005\|LIB10006\|LIB10007\|LIB10008\|LIB10009\|LIB10010\|LIB10011\|LIB10016\|LIB10017\|LIB10018\|LIB10019\|LIB9980\|LIB9981\|LIB9982\|LIB9983\|LIB9996\|LIB9997\|LIB9998\|LIB9999'
# Plate C:
#plateToDo='LIB11214\|LIB11215\|LIB11216\|LIB11217\|LIB11218\|LIB11219\|LIB11220\|LIB11221\|LIB11222\|LIB11223\|LIB11224\|LIB11225\|LIB11226\|LIB11228\|LIB11229\|LIB11230\|LIB11231\|LIB11232\|LIB11233\|LIB11234\|LIB11235\|LIB11236\|LIB11238\|LIB11239\|LIB11240\|LIB11241\|LIB11243\|LIB11244'
#plateToDo='LIB11178\|LIB11179\|LIB11180\|LIB11181\|LIB11186\|LIB11187\|LIB11188\|LIB11189\|LIB11190\|LIB11191\|LIB11192\|LIB11193\|LIB11194\|LIB11195\|LIB11196\|LIB11197\|LIB11198\|LIB11199\|LIB11200\|LIB11202\|LIB11203\|LIB11204\|LIB11206\|LIB11207\|LIB11208\|LIB11209\|LIB11210\|LIB11211\|LIB11212'
#plateToDo='LIB11246\|LIB11247\|LIB11248\|LIB11249\|LIB11250\|LIB11252\|LIB11253\|LIB11254\|LIB11257\|LIB11258\|LIB11259\|LIB11260\|LIB11261\|LIB11262\|LIB11263\|LIB11264\|LIB11265\|LIB11266\|LIB11267\|LIB11269\|LIB11270\|LIB11271\|LIB11272\|LIB11273'
# PlateD:
#plateToDo='LIB10410\|LIB10403\|LIB10404\|LIB10409\|LIB10406\|LIB10407\|LIB10495\|LIB10411\|LIB10412\|LIB10413\|LIB10414\|LIB10415\|LIB10416\|LIB10417\|LIB10418\|LIB10419\|LIB10420\|LIB10421\|LIB10422\|LIB10427\|LIB10428\|LIB10429\|LIB10430\|LIB10431\|LIB10432\|LIB10433\|LIB10434\|LIB10435\|LIB10436\|LIB10437\|LIB10438'
#plateToDo='LIB10439\|LIB10440\|LIB10441\|LIB10442\|LIB10443\|LIB10444\|LIB10445\|LIB10446\|LIB10447\|LIB10448\|LIB10449\|LIB10450\|LIB10399\|LIB10451\|LIB10453\|LIB10455\|LIB10456\|LIB10457\|LIB10458\|LIB10459\|LIB10461\|LIB10462\|LIB10463\|LIB10464\|LIB10465\|LIB10466'
#plateToDo='LIB10471\|LIB10472\|LIB10473\|LIB10474\|LIB10475\|LIB10476\|LIB10477\|LIB10478\|LIB10479\|LIB10480\|LIB10481\|LIB10482\|LIB10487\|LIB10488\|LIB10489\|LIB10490\|LIB10491\|LIB10492\|LIB10493'
# PlateE:
#plateToDo='LIB10289\|LIB10290\|LIB10291\|LIB10292\|LIB10217\|LIB10218\|LIB10219\|LIB10220\|LIB10226\|LIB10227\|LIB10228\|LIB10225\|LIB10230\|LIB10232\|LIB10237\|LIB10313\|LIB10233\|LIB10234\|LIB10244\|LIB10248\|LIB10240\|LIB10243\|LIB10314\|LIB10315\|LIB10242\|LIB10261\|LIB10263\|LIB10264\|LIB10247\|LIB10265\|LIB10266\|LIB10267'
#plateToDo='LIB10249\|LIB10250\|LIB10251\|LIB10252\|LIB10254\|LIB10255\|LIB10256\|LIB10257\|LIB10258\|LIB10259\|LIB10260\|LIB10269\|LIB10270\|LIB10271\|LIB10272\|LIB10273\|LIB10274\|LIB10275\|LIB10276\|LIB10278\|LIB10279\|LIB10280\|LIB10281\|LIB10282\|LIB10283\|LIB10284\|LIB10317\|LIB10318'
#plateToDo='LIB10285\|LIB10286\|LIB10287\|LIB10288\|LIB10293\|LIB10294\|LIB10295\|LIB10296\|LIB10297\|LIB10298\|LIB10299\|LIB10300\|LIB10301\|LIB10302\|LIB10303\|LIB10304\|LIB10305\|LIB10306\|LIB10307\|LIB10308\|LIB10309\|LIB10310\|LIB10311\|LIB10312'
# PlateF:
#plateToDo='LIB10912\|LIB10913\|LIB10914\|LIB10915\|LIB10916\|LIB10917\|LIB10918\|LIB10919\|LIB10920\|LIB10921\|LIB10922\|LIB10927\|LIB10928\|LIB10929\|LIB10930\|LIB10931\|LIB10932\|LIB10933\|LIB10934\|LIB10935\|LIB10936\|LIB10937\|LIB10938\|LIB11007'
#plateToDo='LIB10947\|LIB10948\|LIB10949\|LIB10950\|LIB10951\|LIB10952\|LIB10953\|LIB10955\|LIB10956\|LIB10957\|LIB10958\|LIB10959\|LIB10960\|LIB10961\|LIB10962\|LIB10963\|LIB10964\|LIB10965\|LIB10966\|LIB10967\|LIB10968\|LIB10969\|LIB10970'
#plateToDo='LIB10975\|LIB10976\|LIB10977\|LIB10978\|LIB10979\|LIB10980\|LIB10981\|LIB10982\|LIB10983\|LIB10984\|LIB10985\|LIB10986\|LIB10987\|LIB10988\|LIB10990\|LIB11008\|LIB10991\|LIB10992\|LIB10993\|LIB10995\|LIB10996\|LIB10997\|LIB10998\|LIB10999\|LIB11000\|LIB11009'
#plateToDo='LIB10954'
# PlateG:
#plateToDo='LIB11394\|LIB11395\|LIB11396\|LIB11397\|LIB11398\|LIB11399\|LIB11400\|LIB11401\|LIB11403\|LIB11404\|LIB11405\|LIB11402\|LIB11406\|LIB11407\|LIB11408\|LIB11409\|LIB11410\|LIB11411\|LIB11412\|LIB11413\|LIB11414\|LIB11415\|LIB11417\|LIB11490\|LIB11418\|LIB11419\|LIB11420\|LIB11421\|LIB11422\|LIB11423\|LIB11425\|LIB11491'
#plateToDo='LIB11458\|LIB11459\|LIB11460\|LIB11461\|LIB11462\|LIB11463\|LIB11464\|LIB11465\|LIB11466\|LIB11467\|LIB11468\|LIB11469\|LIB11470\|LIB11473\|LIB11494\|LIB11474\|LIB11475\|LIB11476\|LIB11477\|LIB11478\|LIB11479\|LIB11481\|LIB11495\|LIB11483\|LIB11485'
#plateToDo='LIB11426\|LIB11427\|LIB11428\|LIB11429\|LIB11430\|LIB11431\|LIB11432\|LIB11433\|LIB11434\|LIB11435\|LIB11436\|LIB11437\|LIB11438\|LIB11439\|LIB11440\|LIB11441\|LIB11442\|LIB11443\|LIB11444\|LIB11445\|LIB11446\|LIB11447\|LIB11449\|LIB11492\|LIB11450\|LIB11451\|LIB11452\|LIB11453\|LIB11454\|LIB11455\|LIB11493'
# PlateH:|
#plateToDo='LIB15585\|LIB15586\|LIB15587\|LIB15589\|LIB15590\|LIB15592\|LIB15593\|LIB15594\|LIB15596\|LIB15599\|LIB15601\|LIB15602\|LIB15603\|LIB15604\|LIB15605\|LIB15606\LIB15607\|LIB15608\|LIB15609\|LIB15610\|LIB15611\|LIB15612\|LIB15613\|LIB15614\|LIB15616'
#plateToDo='LIB15617\|LIB15618\|LIB15619\|LIB15620\|LIB15621\|LIB15622\|LIB15623\|LIB15624\|LIB15625\|LIB15626\|LIB15627\|LIB15628\|LIB15629\|LIB15630\|LIB15631\|LIB15632\|LIB15633\|LIB15634\|LIB15635\|LIB15636\|LIB15637\|LIB15638\|LIB15641\|LIB15642\|LIB15643\|LIB15644\|LIB15645\|LIB15646\|LIB15647\|LIB15648'
#plateToDo='LIB15674'
# PlateH + 12 Plate20 samples:
#plateToDo='LIB15657\|LIB15660\|LIB15661\|LIB15662\|LIB15665\|LIB15666\|LIB15667\|LIB15670\|LIB15671\|LIB15672\|LIB15675\|LIB15676\|LIB15677\|LIB15678\|LIB15679\|LIB15680' # PlateH
#plateToDo='LIB17351\|LIB17352\|LIB17353\|LIB17354\|LIB17355\|LIB17356\|LIB17357\|LIB17358\|LIB17359\|LIB17360\|LIB17361\|LIB17362'	# Plate20 (these are in PlateN dir)
# PlateI:
#plateToDo='LIB15289\|LIB15290\|LIB15291\|LIB15292\|LIB15293\|LIB15294\|LIB15296\|LIB15297\|LIB15298\|LIB15299\|LIB15300\|LIB15301\|LIB15302\|LIB15303\|LIB15304\|LIB15305\|LIB15306\|LIB15307\|LIB15308\|LIB15309\|LIB15310\|LIB15311\|LIB15312\|LIB15313\|LIB15314\|LIB15315\|LIB15316\|LIB15317\|LIB15318\|LIB15319\|LIB15320'
#plateToDo='LIB15321\|LIB15322\|LIB15323\|LIB15324\|LIB15325\|LIB15327\|LIB15328\|LIB15329\|LIB15330\|LIB15331\|LIB15332\|LIB15333\|LIB15334\|LIB15335\|LIB15336\|LIB15337\|LIB15338\|LIB15339\|LIB15340\|LIB15341\|LIB15342\|LIB15343\|LIB15344\|LIB15345\|LIB15346\|LIB15347\|LIB15348\|LIB15349\|LIB15350\|LIB15351\|LIB15352'
#plateToDo='LIB15353\|LIB15354\|LIB15355\|LIB15356\|LIB15357\|LIB15358\|LIB15359\|LIB15361\|LIB15362\|LIB15363\|LIB15364\|LIB15365\|LIB15366\|LIB15367\|LIB15368\|LIB15369\|LIB15370\|LIB15371\|LIB15372\|LIB15373\|LIB15374\|LIB15375\|LIB15376\|LIB15377\|LIB15378\|LIB15379\|LIB15380\|LIB15381\|LIB15382\|LIB15383\|LIB15384'
# PlateJ:
#plateToDo='LIB16072\|LIB16073\|LIB16074\|LIB16075\|LIB16076\|LIB16077\|LIB16078\|LIB16079\|LIB16080\|LIB16081\|LIB16082\|LIB16083\|LIB16084\|LIB16085\|LIB16087\|LIB16088\|LIB16090\|LIB16091\|LIB16092\|LIB16093\|LIB16094\|LIB16095\|LIB16096\|LIB16097\|LIB16098\|LIB16100\|LIB16101\|LIB16102\|LIB16103'
#plateToDo='LIB16104\|LIB16105\|LIB16106\|LIB16107\|LIB16108\|LIB16109\|LIB16111\|LIB16112\|LIB16114\|LIB16115\|LIB16116\|LIB16117\|LIB16118\|LIB16119\|LIB16120\|LIB16121\|LIB16122\|LIB16123\|LIB16124\|LIB16125\|LIB16126\|LIB16127\|LIB16128\|LIB16129\|LIB16130\|LIB16131\|LIB16132\|LIB16133\|LIB16135'
#plateToDo='LIB16136\|LIB16138\|LIB16139\|LIB16140\|LIB16141\|LIB16142\|LIB16143\|LIB16144\|LIB16145\|LIB16146\|LIB16147\|LIB16148\|LIB16149\|LIB16150\|LIB16151\|LIB16152\|LIB16153\|LIB16154\|LIB16155\|LIB16156\|LIB16157\|LIB16158\|LIB16159\|LIB16160\|LIB16162\|LIB16163\|LIB16164\|LIB16165\|LIB16166\|LIB16167'
#plateToDo='LIB16165'
# PlateK:
#plateToDo='LIB15880\|LIB15881\|LIB15882\|LIB15884\|LIB15885\|LIB15886\|LIB15887\|LIB15888\|LIB15889\|LIB15890\|LIB15892\|LIB15893\|LIB15894\|LIB15903\|LIB15953\|LIB15954\|LIB15956\|LIB15957\|LIB15958\|LIB15959\|LIB15929\|LIB15945\|LIB15946\|LIB15947\|LIB15948\|LIB15950\|LIB15951'
#plateToDo='LIB15936\|LIB15937\|LIB15938\|LIB15939\|LIB15940\|LIB15941\|LIB15942\|LIB15960\|LIB15961\|LIB15963\|LIB15964\|LIB15965\|LIB15966\|LIB15967\|LIB15968\|LIB15969\|LIB15970\|LIB15971\|LIB15972\|LIB15973\|LIB15974'
#plateToDo='LIB15904\|LIB15905\|LIB15906\|LIB15908\|LIB15909\|LIB15910\|LIB15911\|LIB15912\|LIB15913\|LIB15914\|LIB15915\|LIB15916\|LIB15918\|LIB15919\|LIB15920\|LIB15921\|LIB15923\|LIB15924\|LIB15928\|LIB15930\|LIB15932\|LIB15933\|LIB15934\|LIB15935'
# PlateL:
#plateToDo='LIB16168\|LIB16170\|LIB16171\|LIB16172\|LIB16173\|LIB16174\|LIB16175\|LIB16169\|LIB16257\|LIB16258\|LIB16259\|LIB16260\|LIB16263\|LIB16176\|LIB16177\|LIB16178\|LIB16179\|LIB16180\|LIB16182\|LIB16183\|LIB16200\|LIB16201\|LIB16202\|LIB16203\|LIB16204\|LIB16205\|LIB16206\|LIB16207'
#plateToDo='LIB16208\|LIB16209\|LIB16210\|LIB16212\|LIB16213\|LIB16214\|LIB16215\|LIB16216\|LIB16217\|LIB16218\|LIB16219\|LIB16220\|LIB16221\|LIB16222\|LIB16223\|LIB16248\|LIB16249\|LIB16250\|LIB16251\|LIB16252\|LIB16253\|LIB16254\|LIB16255\|LIB16224\|LIB16225\|LIB16226\|LIB16227\|LIB16228\|LIB16229\|LIB16230\|LIB16231'
# PlateL + 5 Plate20 samples:
#plateToDo='LIB16232\|LIB16233\|LIB16234\|LIB16235\|LIB16236\|LIB16237\|LIB16238\|LIB16239\|LIB16240\|LIB16241\|LIB16242\|LIB16243\|LIB16244\|LIB16245\|LIB16246\|LIB16247'	# PlateL
#plateToDo='LIB15419\|LIB15424\|LIB15441\|LIB15443\|LIB15445'	# Plate20 (these are in PlateN dir)
# PlateM:
#plateToDo='LIB16530\|LIB16531\|LIB16532\|LIB16533\|LIB16534\|LIB16535\|LIB16536\|LIB16537\|LIB16538\|LIB16539\|LIB16540\|LIB16541\|LIB16542\|LIB16543\|LIB16544\|LIB16545\|LIB16546\|LIB16547\|LIB16548\|LIB16549\|LIB16550\|LIB16551\|LIB16552\|LIB16554\|LIB16555\|LIB16556\|LIB16557\|LIB16558\|LIB16559\|LIB16560\|LIB16561'
#plateToDo='LIB16562\|LIB16563\|LIB16564\|LIB16565\|LIB16566\|LIB16567\|LIB16568\|LIB16570\|LIB16571\|LIB16572\|LIB16573\|LIB16574\|LIB16575\|LIB16576\|LIB16578\|LIB16579\|LIB16580\|LIB16581\|LIB16582\|LIB16583\|LIB16584\|LIB16585\|LIB16586\|LIB16587\|LIB16588\|LIB16589\|LIB16590\|LIB16591\|LIB16592\|LIB16593'
#plateToDo='LIB16594\|LIB16595\|LIB16596\|LIB16597\|LIB16598\|LIB16599\|LIB16600\|LIB16602\|LIB16603\|LIB16604\|LIB16605\|LIB16606\|LIB16607\|LIB16608\|LIB16609\|LIB16610\|LIB16611\|LIB16612\|LIB16613\|LIB16614\|LIB16615\|LIB16616\|LIB16617\|LIB16618\|LIB16619\|LIB16620\|LIB16621\|LIB16622\|LIB16623\|LIB16624'
plateToDo='LIB16572'
# PlateN:
#plateToDo='LIB15449\|LIB15451\|LIB15452\|LIB15453\|LIB15454\|LIB15455\|LIB15456\|LIB15457\|LIB15458\|LIB15459\|LIB15460\|LIB15461\|LIB15463\|LIB15465\|LIB15467\|LIB15469\|LIB15470\|LIB15472\|LIB15473\|LIB15475\|LIB15476\|LIB15477\|LIB15478\|LIB15479\|LIB15480\|LIB15393\|LIB15395\|LIB15397\|LIB15398\|LIB15400'
#plateToDo='LIB15385\|LIB15386\|LIB15387\|LIB15388\|LIB15389\|LIB15391\|LIB15401\|LIB15403\|LIB15404\|LIB15405\|LIB15409\|LIB15410\|LIB15411\|LIB15412\|LIB15413\|LIB15415\|LIB15425\|LIB15427\|LIB15429\|LIB15431\|LIB15432\|LIB15433\|LIB15434\|LIB15436\|LIB15439'

# Poor captures:
# Plates A --> G:
#bamFileLoctn=/tgac/workarea/group-cg/baileyp/WheatLoLa/35.ExonCapture_myIWGSC_v2_chrU_ref/Plate[ACEFG]	# Plates [ACEFG]:
#bamFileLoctn=/tgac/scratch/baileyp/WheatLoLa/35.ExonCapture_myIWGSC_v2_chrU_ref/Plate[BD]	# Plates [BD]:
# List of captures (On-targReadPairs < 5m) - adding these to a separate run of the 19 5m -10m samples - T = 27:
#plateToDo='LIB10612\|LIB10614\|LIB10615\|LIB10637\|LIB10672\|LIB10674\|LIB10675\|LIB9956\|LIB9957\|LIB9958\|LIB9959\|LIB10454\|LIB10460\|LIB10231\|LIB10235\|LIB10241\|LIB11002\|LIB11484\|LIB11486\|LIB10943\|LIB10944\|LIB10945\|LIB10946\|LIB10701'
#plateToDo='LIB9944\|LIB9946\|LIB9947\|LIB9956\|LIB9957\|LIB9958\|LIB9959\|LIB10454\|LIB10460'
#plateToDo='LIB9956\|LIB9957\|LIB9958\|LIB9959\|LIB10454\|LIB10460'
#plateToDo='LIB10607'
#List of captures (On-targReadPairs ~5m - ~10m):
#plateToDo='LIB10612\|LIB10614\|LIB10615\|LIB10637\|LIB10672\|LIB10674\|LIB10675\|LIB10231\|LIB10235\|LIB10241\|LIB11002\|LIB11484\|LIB11486\|LIB10943\|LIB10944\|LIB10945\|LIB10946\|LIB10701'
#plateToDo='LIB9956\|LIB9957\|LIB9958\|LIB9959\|LIB10454\|LIB10460'

#List of captures (On-targReadPairs < ~10m - ~15m):
#plateToDo='LIB10673\|LIB10674\|LIB10675\|LIB10676\|LIB10677\|LIB10679\|LIB10689\|LIB10695\|LIB10628\|LIB10687\|LIB10994\|LIB11457\|LIB11201\|LIB11227\|LIB11251\|LIB11255\|LIB11256\|LIB11268\|LIB11205\|LIB11213\|LIB11237\|LIB11242\|LIB11245'
#plateToDo='LIB9940\|LIB9941\|LIB9942\|LIB9943\|LIB10012\|LIB10013\|LIB10014\|LIB10015'

#List of captures (On-targReadPairs < ~15m - ~20m):
#bamFileLoctn=/tgac/workarea/group-cg/baileyp/WheatLoLa/35.ExonCapture_myIWGSC_v2_chrU_ref/Plate[ACEFG]	# Plates [ACEFG]:
#plateToDo='LIB10664\|LIB10665\|LIB10666\|LIB10667|LIB10668\|LIB10669\|LIB10671\|LIB10644\|LIB10645\|LIB10646\|LIB10647\|LIB10648\|LIB10649\|LIB10650\|LIB10651\|LIB10673\|LIB10943\|LIB10944\|LIB10945\|LIB10946\|LIB10701\|LIB10702\|LIB10705\|LIB10706'
#plateToDo='LIB10667\|LIB10668'

# Poor captures - plates H to N:
#bamFileLoctn=/tgac/workarea/group-cg/baileyp/WheatLoLa/39.ExonCapture_PlateH2N/Plate[IJKLM]	# Plates [IJKLM]:
#bamFileLoctn=/tgac/scratch/baileyp/WheatLoLa/39.ExonCapture_PlateH2N/Plate[HN]	# Plates [HN]:
#List of captures (On-targReadPairs < 5 million) - adding these to a separate run of the 19 5-10m samples - T = 27:
#plateToDo='LIB15925\|LIB15926\|LIB15927\|LIB16227\|LIB16099\|LIB16148\|LIB16150\|LIB16158\|LIB15891\|LIB15895\|LIB16220\|LIB16569'	# Plate IJKLM
#plateToDo='LIB15597\|LIB15658\|LIB15437\|LIB15640\|LIB15659\|LIB15664\|LIB15406\|LIB15407\|LIB15408\|LIB15428\|LIB15430\|LIB15435\|LIB15478\|LIB15417\|LIB15448'	# Plate H+N
#plateToDo='LIB15435'

#List of captures (On-targReadPairs  ~5m - ~10m):  
#plateToDo='LIB16099\|LIB16148\|LIB16150\|LIB16158\|LIB15891\|LIB15895\|LIB16220\|LIB16569' # PlateIJKLM
#plateToDo='LIB15640\|LIB15659\|LIB15664\|LIB15406\|LIB15407\|LIB15408\|LIB15428\|LIB15430\|LIB15478\|LIB15417\|LIB15448' # PlateH+N

for chrArm in ${chrArms[@]}; do

	echo $chrArm
	
	# Only make directory for the first run: 
	if [[ ! -d $chrArm ]]; then mkdir $chrArm; fi
	
	cd $chrArm

	bedfile=${splitRefBedfileLocatn}/${chrArm}.fa.bed
	echo $bedfile
	

	tail -n +2 $fileList | grep $plateToDo |
	while read line; do

		pathsToLibs=`echo $line | cut -d ' ' -f 5`

		sampleDir=`basename $pathsToLibs`
		#echo $sampleDir
		# NB - newer seqing runs and newer libraries have longer ids
		### 29.6.0215 - need to check this again for the old and new libraries -  seem to need single quotes round the sed cmd
		### 21.7.2015 - could also test this line again for PlatesH-->N - improved the search syntax
		### 14.9.2015 - some old test capture samples only have 3 digits after "Sample_":
		samplePrefix=`echo $sampleDir | sed 's/Sample_[0-9]\{3,\}_//' | sed 's/_LDI[0-9]\{4,\}//'`
		# For PlatesH-N, e.g. LIB15321_LDI13200
		#samplePrefix=`echo $sampleDir | sed 's/_LDI[0-9]\{4,\}//'`
		echo $samplePrefix
		#echo $bamFileName
		echo $bamFileLoctn/${sampleDir}/$bamFileName

		bsub -J ${samplePrefix}_$chrArm -n 1 -q normal -R "rusage [mem=7500]" -oo ${chrArm}_${samplePrefix}_intersect_samtools_index.log "
		# NB - ||elized mpileup doesn't work with these bam files, with or without -header flag
		source bedtools-2.17.0; bedtools intersect -header \
		-abam $bamFileLoctn/${sampleDir}/$bamFileName  \
		-b $bedfile \
		> ${samplePrefix}_sample.bam
		# Index the bam file (only required for the ||elized version of run-mpileup.py):
		#source /tgac/software/production/bin/samtools-0.1.18
		#samtools index ${samplePrefix}_sample.bam
		"
	done
	
	cd ../
done
***string_marker_for_commenting_out_code***


# An easy way of switching jobs to different queues:- e.g.:
#queues=(TempProject1 Test128 Prod256 Test256 Prod128)
#while true; do for queue in ${queues[@]}; do echo " "; echo " "; echo " "; echo " "; echo " "; echo " "; echo " "; echo " "; sleep 5; echo $queue; for jobId in {8631857..8632522}; do bmod -q $queue $jobId; done; done; done








#############################
#Set parameters for run here:
#############################
runMpileupVersion=~baileyp/ProgramFiles/maps_fortgac/run-mpileup.py										# Other version: run-mpileup[_-A].py
#runMpileupVersion=/usr/users/ga002/baileyp/ProgramFiles/mpileup-tools-master/beta-run-mpileup.py		# ||elized version of run_mpileup - beta-run-mpileup[_-A].py



#MAPS-part1:
jobId=maps				
minLib=20				#***************NEED TO SET THIS EVERY TIME***************
minCov=10
hetOneMinPer=20			# default = 20
maxCov=10000			# default = 2000

#MAPS-part2:
hetMinPer=10			# Running with -p=20 (original value I used), 15, 25 and 40
#homMinCov=3

:<<***string_marker_for_commenting_out_code***

# Remember:
# Start MAPS run in the UV2 /scratch/ dir and save subsequent outputs to there - writing to /scratch is much quicker 

cpu=4				# Set to multiples of 4 to fit into a UV2 nodes that have multiples of 4 cpus; used 48 on the UV2; 24 on the Cluster; also trying 96 
mem_mb=6500			# Required for bsub - specify 210000 mb(actually ~203 GB) for 32 samples; NB - it seems like 3.6 GB per thread (e.g. 700GB mem for 192 threads); for /dev/shm, use 3500000; 6500 in chunk mode - NB - TGAC_cadenza_U looks to need up to 12,000 
mem_gb=6.5			# Required for qsub - specify 210gb (actually ~203 GB) for 32 samples; for /dev/shm, use 3500
dev_shm=			# e.g. /dev/shm/PlateE/MAPS_bams_to_use_LIB10289_LIB10267/ (note the trailing "/" char and any dirs need to exist!) ; puts files and everything into memory; leave blank for a normal run


# Running MAPS separately for each chr arm:
for chrArm in ${chrArms[@]}; do

	echo $chrArm
	cd $chrArm
	
	chrArmReferenceFile=${splitRefBedfileLocatn}/${chrArm}.fa
	echo $chrArmReferenceFile

#:<<***string_marker_for_commenting_out_code***
	
	### NBNB - can't get the nested qsubs to work cleanly! 
	### 28.4.2015 - trialing use of the ampersand instead after the maps-paret2 step; also adding "sleep 10" just afterwards - 
	### this might allow error messages to be collected by the logfile (?)
	#echo "cd $PWD;
	bsub -J maps_$chrArm -q Prod128 -n $cpu -R "rusage[mem=${mem_mb}]" -oo maps_run_-l${minLib}.log "

	# For the uv's, need to source hpcore:
	source /tgac/software/testing/bin/hpccore-5
	source /tgac/software/testing/bin/python-2.7.5
	# ||elized version:
	#python $runMpileupVersion \
	#	-t $cpu \
	#	-n _sample.bam \
	#	-d 8000 \
	#	-r $chrArmReferenceFile \
	#	-o ${dev_shm}run_mpileup.out \
	#	-s /tgac/software/production/samtools/0.1.18/x86_64/bin/samtools
	# Non-||elized version:
#:<<***string_marker_for_commenting_out_code***
	python $runMpileupVersion \
		-d 8000 \
		-r $chrArmReferenceFile \
		-o run_mpileup.out \
		-s /tgac/software/production/samtools/0.1.18/x86_64/bin/samtools
	python ~baileyp/ProgramFiles/mpileup-tools-master/beta-mpileup-parser.py \
		-t $cpu \
		-f ${dev_shm}run_mpileup.out \
		-o ${dev_shm}parsed_run_mpileup.out
	python ~baileyp/ProgramFiles/maps_fortgac/maps-part1-v2.py \
		-t $cpu \
		-f ${dev_shm}parsed_run_mpileup.out \
		-o ${dev_shm}maps-part1_l${minLib}_v${minCov}_i${hetOneMinPer}.txt \
		-m m \
		-l $minLib \
		-v $minCov \
		-i $hetOneMinPer \
		--maxCov $maxCov
#***string_marker_for_commenting_out_code***
	hetMinCovs=(1 2 3 4 5 6 7 8 9 10 11 12)
	homMinCov=2		# was 3
### 3.2.2016 - In future just loop through several values of homMinCov as well: 1,2,3,4,5,6(,7,8,9,10)
### Actually to get the homMinCov=1 to work here, need to backslash it.    
	
	for hetMinCov in \${hetMinCovs[@]}; do

		echo \"hetMinCov=\$hetMinCov	hetMinPer=$hetMinPer	homMinCov=\$homMinCov\"
	
		# bsubing these to run in ||el (RUN TIME: 3.5 hours):
		bsub -J maps2_c\${hetMinCov} -q Prod128 -n 1 -R "rusage[mem=1000]" -oo maps-part2_-l${minLib}_-c${minCov}_-d\${hetMinCov}_-p${hetMinPer}_-s\$homMinCov.log \"

		python ~baileyp/ProgramFiles/maps_fortgac/maps-part2-v2.py \
		-f ${dev_shm}maps-part1_l${minLib}_v${minCov}_i${hetOneMinPer}.txt \
		-o maps-part2_-l${minLib}_-c${minCov}_-d\${hetMinCov}_-p${hetMinPer}_-s\$homMinCov.txt \
		-m m \
		-l ${minLib} \
		-c $minCov \
		-d \$hetMinCov \
		-p $hetMinPer \
		-s \$homMinCov
		\"
	done
	homMinCov=3		# Access this variable like so: \$homMinCov!

	for hetMinCov in \${hetMinCovs[@]}; do

		echo \"hetMinCov=\$hetMinCov	hetMinPer=$hetMinPer	homMinCov=\$homMinCov\"

		bsub -J maps2_c\${hetMinCov} -q Prod128 -n 1 -R "rusage[mem=1000]" -oo maps-part2_-l${minLib}_-c${minCov}_-d\${hetMinCov}_-p${hetMinPer}_-s\$homMinCov.log \"

		python ~baileyp/ProgramFiles/maps_fortgac/maps-part2-v2.py \
		-f ${dev_shm}maps-part1_l${minLib}_v${minCov}_i${hetOneMinPer}.txt \
		-o maps-part2_-l${minLib}_-c${minCov}_-d\${hetMinCov}_-p${hetMinPer}_-s\$homMinCov.txt \
		-m m \
		-l ${minLib} \
		-c $minCov \
		-d \$hetMinCov \
		-p $hetMinPer \
		-s \$homMinCov
		\"
	done
	#if [[ -a run_mpileup.out ]]; then rm run_mpileup.out; fi
	#if [[ -a parsed_run_mpileup.out ]]; then rm parsed_run_mpileup.out; fi
	"
	#| qsub -N $jobId -q Prod -l select=1:ncpus=$cpu:mem=${mem_gb}gb -j eo
#***string_marker_for_commenting_out_code***
	
	cd ../

done
***string_marker_for_commenting_out_code***





# Concatenate maps-part2_*.txt files from all MAPS runs for each hetMinCov value (RUN TIME: 8 minutes):
# Run command: bsub -oo maps-part2_-lX_-c10_-dall_-p10_-sall_platesABCD.log "35.runMAPS_InChrArmChunks.sh"
mm=			# NB - for the multiply mapped MAPS runs prepend the "MAPS_bams_to_use*" directory names with "mm" - leave blank for the main MAPS runs 	
#:<<***string_marker_for_commenting_out_code***
hetMinCovs=(1 2 3 4 5 6 7 8 9 10 11 12)		# (1 2 3 4 5 6 7 8 9 10 11 12)
homMinCovs=(2)							# (3 4)
for homMinCov in ${homMinCovs[@]}; do
	for hetMinCov in ${hetMinCovs[@]}; do

		echo "hetMinCov=$hetMinCov	homMinCov=$homMinCov"
		# Count to check the number of files in each set e.g. :
		#ls /tgac/workarea/group-cg/baileyp/WheatLoLa/39.ExonCapture_PlateH2N/Plate[JKL]/MAPS_bams_to_use_*/*/maps-part2_-l*_-c10_-d${hetMinCov}_-p10_-s${homMinCov}.txt | grep -v '^hetMinCov' | wc -l

		# Prepare a concatenated file without any header lines.
		# Also removing the filenames from the tail command.
		# Prepare final MAPS output for plates A --> N (NB - there is now a file for each chr arm)
		bsub -J mps_cat${hetMinCov} -q Prod128 -n 1 -oo ${mm}maps-part2_-lX_-c10_-d${hetMinCov}_-p10_-s${homMinCov}_TABCDEFGHIJKLMN_poorA-G_poorH-N_cat.log "

		# Test captures (the libraries indicated are present in other parts of the data so are removed here):
### NBNB - APRIL 2016 - THESE GREPS WOULD NEED TO BE CHECKED TO REMOVE THE 6 DUPLICATE SAMPELS AND LIB5773 CAN BE KEPT IN]
		tail -n +2 /tgac/workarea/group-cg/baileyp/WheatLoLa/35.ExonCapture_myIWGSC_v2_chrU_ref/54_test_samples/${mm}MAPS_bams_to_use*/*/maps-part2_-l*_-c10_-d${hetMinCov}_-p10_-s${homMinCov}.txt \
		| grep -v 'LIB577[2345679]\|LIB5780\|LIB5789' \
		| grep -v '==> ' \
		> maps-part2_-lX_-c10_-d${hetMinCov}_-p10_-s${homMinCov}_TABCDEFGHIJKLMN_poorA-G_poorH-N.txt

		# Plates A to G:
		tail -n +2 /tgac/workarea/group-cg/baileyp/WheatLoLa/35.ExonCapture_myIWGSC_v2_chrU_ref/Plate[ACEFG]/${mm}MAPS_*/*/maps-part2_-l*_-c10_-d${hetMinCov}_-p10_-s${homMinCov}.txt \
		| grep -v '==> ' \
		>> maps-part2_-lX_-c10_-d${hetMinCov}_-p10_-s${homMinCov}_TABCDEFGHIJKLMN_poorA-G_poorH-N.txt
		tail -n +2 /tgac/scratch/baileyp/WheatLoLa/35.ExonCapture_myIWGSC_v2_chrU_ref/Plate[BD]/${mm}MAPS_*/*/maps-part2_-l*_-c10_-d${hetMinCov}_-p10_-s${homMinCov}.txt \
		| grep -v '==> ' \
		>> maps-part2_-lX_-c10_-d${hetMinCov}_-p10_-s${homMinCov}_TABCDEFGHIJKLMN_poorA-G_poorH-N.txt
 
		# Plates H --> N:			
		tail -n +2 /tgac/workarea/group-cg/baileyp/WheatLoLa/39.ExonCapture_PlateH2N/Plate[JKLM]/${mm}MAPS_*/*/maps-part2_-l*_-c10_-d${hetMinCov}_-p10_-s${homMinCov}.txt \
		| grep -v '==> ' \
		>> maps-part2_-lX_-c10_-d${hetMinCov}_-p10_-s${homMinCov}_TABCDEFGHIJKLMN_poorA-G_poorH-N.txt
		tail -n +2 /tgac/scratch/baileyp/WheatLoLa/39.ExonCapture_PlateH2N/Plate[HIN]/${mm}MAPS_*/*/maps-part2_-l*_-c10_-d${hetMinCov}_-p10_-s${homMinCov}.txt \
		| grep -v 'LIB15478|\LIB16099'
		| grep -v '==> ' \
		>> maps-part2_-lX_-c10_-d${hetMinCov}_-p10_-s${homMinCov}_TABCDEFGHIJKLMN_poorA-G_poorH-N.txt 
		
		# Poor captures for plates A to G:
		# 0-5k: prefiltered to keep only LIB9944\|LIB9946\|LIB9947:
		tail -n +2 /tgac/workarea/group-cg/baileyp/WheatLoLa/35.ExonCapture_myIWGSC_v2_chrU_ref/PlateA-G_poor_captures_MAPS/${mm}less5m_On-targReadPairs/*/maps-part2_-l*_-c10_-d${hetMinCov}_-p10_-s${homMinCov}.txt \
		| grep -v 'LIB10612\|LIB10614\|LIB10615\|LIB10637\|LIB10672\|LIB10674\|LIB10675\|LIB9956\|LIB9957\|LIB9958\|LIB9959\|LIB10454\|LIB10460\|LIB10231\|LIB10235\|LIB10241\|LIB11002\|LIB11484\|LIB11486\|LIB10943\|LIB10944\|LIB10945\|LIB10946\|LIB10701' \
		| grep -v '==> ' \
		>> maps-part2_-lX_-c10_-d${hetMinCov}_-p10_-s${homMinCov}_TABCDEFGHIJKLMN_poorA-G_poorH-N.txt
		# 5 - 10k: removing 3 extra samples present elsewhere in the poor captures:
		tail -n +2 /tgac/workarea/group-cg/baileyp/WheatLoLa/35.ExonCapture_myIWGSC_v2_chrU_ref/PlateA-G_poor_captures_MAPS/${mm}5m-10m_On-targReadPairs/*/maps-part2_-l*_-c10_-d${hetMinCov}_-p10_-s${homMinCov}.txt \
		| grep -v 'LIB10674\|LIB10675\|LIB10946\|LIB10943\|LIB10944\|LIB10945' \
		| grep -v '==> ' \
		>> maps-part2_-lX_-c10_-d${hetMinCov}_-p10_-s${homMinCov}_TABCDEFGHIJKLMN_poorA-G_poorH-N.txt
		# 10 - 15k: removing 1 sample present elsewhere in the poor captures:
		tail -n +2 /tgac/workarea/group-cg/baileyp/WheatLoLa/35.ExonCapture_myIWGSC_v2_chrU_ref/PlateA-G_poor_captures_MAPS/${mm}10m-15m_On-targReadPairs/*/maps-part2_-l*_-c10_-d${hetMinCov}_-p10_-s${homMinCov}.txt \
		| grep -v 'LIB10673' \
		| grep -v '==> ' \
		>> maps-part2_-lX_-c10_-d${hetMinCov}_-p10_-s${homMinCov}_TABCDEFGHIJKLMN_poorA-G_poorH-N.txt
		# 15 - 20k
		tail -n +2 /tgac/workarea/group-cg/baileyp/WheatLoLa/35.ExonCapture_myIWGSC_v2_chrU_ref/PlateA-G_poor_captures_MAPS/${mm}15m-20m_On-targReadPairs/*/maps-part2_-l*_-c10_-d${hetMinCov}_-p10_-s${homMinCov}.txt \
		| grep -v '==> ' \
		>> maps-part2_-lX_-c10_-d${hetMinCov}_-p10_-s${homMinCov}_TABCDEFGHIJKLMN_poorA-G_poorH-N.txt

		# Poor captures for plates H to N:				
		# 0-5k: prefiltered to keep only LIB9944\|LIB9946\|LIB9947:
		tail -n +2 /tgac/workarea/group-cg/baileyp/WheatLoLa/39.ExonCapture_PlateH2N/PlateH2N_poor_captures_MAPS/${mm}less_5m_On-targReadPairs/*/maps-part2_-l*_-c10_-d${hetMinCov}_-p10_-s${homMinCov}.txt \
		| grep -v 'LIB15640\|LIB15659\|LIB15664\|LIB15406\|LIB15407\|LIB15408\|LIB15428\|LIB15430\|LIB15478\|LIB15417\|LIB15448\|LIB16099\|LIB16148\|LIB16150\|LIB16158\|LIB15891\|LIB15895\|LIB16220\|LIB16569' \
		| grep -v '==> ' \
		>> maps-part2_-lX_-c10_-d${hetMinCov}_-p10_-s${homMinCov}_TABCDEFGHIJKLMN_poorA-G_poorH-N.txt
		# 5 - 10k: 
		tail -n +2 /tgac/workarea/group-cg/baileyp/WheatLoLa/39.ExonCapture_PlateH2N/PlateH2N_poor_captures_MAPS/${mm}5m-10m_On-targReadPairs/*/maps-part2_-l*_-c10_-d${hetMinCov}_-p10_-s${homMinCov}.txt \
		| grep -v '==> ' \
		>> maps-part2_-lX_-c10_-d${hetMinCov}_-p10_-s${homMinCov}_TABCDEFGHIJKLMN_poorA-G_poorH-N.txt
		"
	done
done
#***string_marker_for_commenting_out_code***
# The maps-part2_*.txt file is used to prepare the Excel file and to get d3, d4 and d5 lines separately (see below)





# Now check that I have the expected number of libraries - just use hetMonCov=6 and homMinCov=3
#awk '{print $7}' maps-part2_-lX_-c10_-d6_-p10_-s3_TABCDEFGHIJKLMN_poorA-G_poorH-N.txt | sort -u | wc -l
# Total number samples = 1244 - correct!





# Now run 35.runMAPS_Pipeline_GetCounts.sh individually on each hetMinCov value/file (RUN TIME: d2 takes 22h others are less except d1 which takes ~12 days!):
:<<***string_marker_for_commenting_out_code***
hetMinCovs=(1 2 3 4 5 6 7 8 9 10 11 12)		# 1 2 3 4 5 6 7 8 9 10 11 12
homMinCovs=(3 4)							# 3 4
for homMinCov in ${homMinCovs[@]}; do
	for hetMinCov in ${hetMinCovs[@]}; do

		printf '%b' "Lane\tLibrary\tDescription\tCadenzaNo.\tTotalSNPs(-d=${hetMinCov})\tHetSNPs\tHomSNPs\tHet/HomRatio\tGA\tCT\tAG\tTC\tAC\tAT\tCA\tCG\tGC\tGT\tTA\tTG\t \
TotalEMS_Ts\tTotalEMS_Ts_%\tTotalNon-EMS_Ts\tTotalNon-EMS_Ts_%\tTotalNon-EMS_Tv\tTotalNon-EMS_Tv_%\t.\t\n \
" >ExonCapture_maps-part2_-d${hetMinCov}_-s${homMinCov}_TABCDEFGHIJKLMN_poorA-G_poorH-N_counts_for_Excel.txt

		# Usage: 35.runMAPS_Pipeline_GetCounts.sh   hetMinCov   homMinCov   MAPS_output   output_for_Excel   sample_info_table
		bsub -J mps_cnts${hetMinCov} -q Prod128 -n 1 -R "rusage[mem=1000]" -oo ExonCapture_maps-part2_-d${hetMinCov}_-s${homMinCov}_TABCDEFGHIJKLMN_poorA-G_poorH-N_counts_for_Excel.log "
		35.runMAPS_Pipeline_GetCounts.sh $hetMinCov $homMinCov \
	    maps-part2_-lX_-c10_-d${hetMinCov}_-p10_-s${homMinCov}_TABCDEFGHIJKLMN_poorA-G_poorH-N.txt 
AMMENDED THIS LINE!!!!!		ExonCapture_maps-part2_-d${hetMinCov}_-s${homMinCov}_TABCDEFGHIJKLMN_poorA-G_poorH-N_counts_for_Excel.txt_just_for_LIB7882 \
		/tgac/workarea/group-cg/baileyp/WheatLoLa/39.ExonCapture_PlateH2N/VEP_analysis_aux_files/ExonCapture_PlateTABCDEFGHIJKLMN.txt
		"
	done
done
***string_marker_for_commenting_out_code***
# NB - after running script, it's important to detect samples with lots of non-EMS mutations and exclude them before the VEP is run.





# Extract hetMinCov and homMinCov values at -d3, -d4 and -d5 etc specifically from the corresponding files with these values for the parameters d3, -d4 and -d5 files.
# These values are needed for:
# 1. the % EMS vs hetMinCov plots
# 2. producing the final data so -d3, -d4, -d5 and -d>=6 data are extracted respectively from -d3 -d4, -d5 and -d>=6 files. 
:<<***string_marker_for_commenting_out_code***
hetMinCovs=(3 4 5)	
homMinCovs=(3)			
for homMinCov in ${homMinCovs[@]}; do
	for hetMinCov in ${hetMinCovs[@]}; do

		# Usage: 35.runMAPS_Pipeline_GetCounts.sh   hetMinCov   homMinCov   MAPS_output   output_for_Excel   sample_info_table
		bsub -J mps_${hetMinCov}_only -q Prod128 -n 1 -R "rusage[mem=1000]" -oo maps-part2_-d${hetMinCov}_-p10_-s${homMinCov}___d${hetMinCov}_s${hetMinCov}_lines_only.log "
	    ls maps-part2_-lX_-c10_-d${hetMinCov}_-p10_-s${homMinCov}_TABCDEFGHIJKLMN_poorA-G_poorH-N.txt
	    cat maps-part2_-lX_-c10_-d${hetMinCov}_-p10_-s${homMinCov}_TABCDEFGHIJKLMN_poorA-G_poorH-N.txt \
	    | awk -v hetMinCov=${hetMinCov} '(\$8 ~ /hom/ && \$10 == hetMinCov) || (\$8 ~ /het/ && \$10 == hetMinCov)' \
	    > maps-part2_-d${hetMinCov}_-p10_-s${homMinCov}___d${hetMinCov}_s${hetMinCov}_lines_only.txt
		"
	done
done
# Extract the remaining SNPs with coverage >= 6 from the -d6 -s3 file:
hetMinCov=(6)
homMinCovs=(3)		
bsub -J mps_more_${hetMinCov} -q Prod128 -n 1 -R "rusage[mem=1000]" -oo maps-part2_-d${hetMinCov}_-p10_-s${homMinCov}___more_d${hetMinCov}__more_s${hetMinCov}_lines.log "
ls maps-part2_-lX_-c10_-d${hetMinCov}_-p10_-s${homMinCov}_TABCDEFGHIJKLMN_poorA-G_poorH-N.txt
cat maps-part2_-lX_-c10_-d${hetMinCov}_-p10_-s${homMinCov}_TABCDEFGHIJKLMN_poorA-G_poorH-N.txt \
| awk '(\$8 ~ /hom/ && \$10 >= 6) || (\$8 ~ /het/ && \$10 >= 6)' \
> maps-part2_-d${hetMinCov}_-p10_-s${homMinCov}___more_d${hetMinCov}__more_s${hetMinCov}_lines.txt
"
***string_marker_for_commenting_out_code***


# Concatenate the above files for coverage of 3, 4, 5, and >= 6:
#cat maps-part2_-d3_-p10_-s3___d3_s3_lines_only.txt \
#maps-part2_-d4_-p10_-s3___d4_s4_lines_only.txt \
#maps-part2_-d5_-p10_-s3___d5_s5_lines_only.txt \
#maps-part2_-d6_-p10_-s3___more_d6__more_s6_lines.txt \
#> maps-part2_3x_4x_5x_more6x_lines.txt
# This file can now be used in SNP annotation





# Paste the outputs for Excel together to make a single table: 
# paste  ExonCapture_maps-part2_-d*_-s3_platesABCDEFG_counts_for_Excel.txt > ExonCapture_maps-part2_-dall_-s3_platesABCDEFG_counts_for_Excel.txt
# 22.9.2015 - waiting for -d1....
#paste  ExonCapture_maps-part2_-d[23456789]_-s3_TABCDEFGHIJKLMN_poorA-G_poorH-N_counts_for_Excel.txt > ExonCapture_maps-part2_-dall_-s3_TABCDEFGHIJKLMN_poorA-G_poorH-N_counts_for_Excel.txt
#paste  ExonCapture_maps-part2_-d[23456789]_-s4_TABCDEFGHIJKLMN_poorA-G_poorH-N_counts_for_Excel.txt > ExonCapture_maps-part2_-dall_-s4_TABCDEFGHIJKLMN_poorA-G_poorH-N_counts_for_Excel.txt 

# Next steps:
# 1. Once data is in Excel, need to remove the ID and info columns for each value of -d after checking that all columns are in sync
# 2. Look at data and prepare to grep out the samples for calculating the stats - low % EMS lines and lines with a high het/hom ratio   