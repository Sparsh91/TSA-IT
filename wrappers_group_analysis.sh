#! /bin/bash

list1='wm_LR_2bk_cor relational_LR_relation relational_LR_match motor_LR_rf motor_LR_lf motor_LR_rh motor_LR_lh language_LR_math language_LR_story gambling_LR_win_event gambling_LR_loss_event gambling_LR_neut_event social_LR_rnd social_LR_mental_resp emotion_LR_neut emotion_LR_fear'

language='story math'
wm='0bk_cor 2bk_cor'
emotion='fear neut'
social='mental_resp rnd'
gambling='neut_event loss_event win_event'
motor='lh rh lf rf'
relational='match relation'

for i in $language;
do
	python group_analysis_part.py language_LR_RL/$i/output_z language_LR_RL/$i/group_analysis
done 

for i in $wm;
do
	python group_analysis_part.py wm_LR_RL/$i/output_z wm_LR_RL/$i/group_analysis
done 
for i in $emotion;
do
	python group_analysis_part.py emotion_LR_RL/$i/output_z emotion_LR_RL/$i/group_analysis
done 
for i in $social;
do
	python group_analysis_part.py social_LR_RL/$i/output_z social_LR_RL/$i/group_analysis
done 
for i in $gambling;
do
	python group_analysis_part.py gambling_LR_RL/$i/output_z gambling_LR_RL/$i/group_analysis
done 
for i in $motor;
do
	python group_analysis_part.py motor_LR_RL/$i/output_z motor_LR_RL/$i/group_analysis
done 
for i in $relational;
do
	python group_analysis_part.py relational_LR_RL/$i/output_z relational_LR_RL/$i/group_analysis
done 
