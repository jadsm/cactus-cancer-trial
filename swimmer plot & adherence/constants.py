google_credentials_path = "/Users/juandelgado/Desktop/Juan/code/creds/gcreds_cactus.json"

# treatment_schedule_classes = {'Targeted_to_immuno':[8,9,1,20,12,7],
#                                         'Immuno_to_targeted':[18],
#                                         'Targeted_only':[3,6,10,5,2],
#                                         'Immuno_only':[15,16,17],
#                                         'double_switch':[4,14,21,11,19,13]}

project_id = 'systemsforecasting'
treatment_schedule_classes = {'T->Io':[8,9,1,20,12,7],
                              'Io->T':[18],
                                'T only':[3,6,10,5,2],
                                'Io only':[15,16,17],
                                'T->Io->T':[4,14,21,11,19,13]}

first_line_treatment_classes = {'Targeted':[8,9,1,20,12,7,3,6,10,5,2,4,14,21,11,19,13],
                              'Io':[18,15,16,17]}

vaf_risk_dictionary = {r'20':'4 - severe risk',#'>20%'
                        r'10': '3 - high risk',
                        r'5': '2 - medium risk',
                        r"3": '1 - lower risk',
                        r'1.5': '0 - not measurable'}