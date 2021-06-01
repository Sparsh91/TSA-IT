#! /bin/bash

list1='wm_LR_2bk_cor relational_LR_relation relational_LR_match motor_LR_rf motor_LR_lf motor_LR_rh motor_LR_lh language_LR_math language_LR_story gambling_LR_win_event gambling_LR_loss_event gambling_LR_neut_event social_LR_rnd social_LR_mental_resp emotion_LR_neut emotion_LR_fear'

language='story math'
wm='0bk_cor 2bk_cor'
emotion='fear neut'
social='mental_resp rnd'
gambling='neut_event loss_event win_event'
motor='lh rh lf rf'
relational='match relation'

mkdir language_LR_RL
for i in $language;
do
	mkdir language_LR_RL/$i
	python wrapper_LR_RL.py language_LR_RL/$i/output language_LR_RL/$i/output_z tfMRI_LANGUAGE $i ev_4dfiles_LR_RL_wrapper
done 

mkdir wm_LR_RL
for i in $wm;
do
	mkdir wm_LR_RL/$i
	python wrapper_LR_RL.py wm_LR_RL/$i/output wm_LR_RL/$i/output_z tfMRI_WM $i ev_4dfiles_LR_RL_wrapper
done 

mkdir emotion_LR_RL
for i in $emotion;
do
	mkdir emotion_LR_RL/$i
	python wrapper_LR_RL.py emotion_LR_RL/$i/output emotion_LR_RL/$i/output_z tfMRI_EMOTION $i ev_4dfiles_LR_RL_wrapper
done 


mkdir social_LR_RL
for i in $social;
do
	mkdir social_LR_RL/$i
	python wrapper_LR_RL.py social_LR_RL/$i/output social_LR_RL/$i/output_z tfMRI_SOCIAL $i ev_4dfiles_LR_RL_wrapper
done 


mkdir gambling_LR_RL
for i in $gambling;
do
	mkdir gambling_LR_RL/$i
	python wrapper_LR_RL.py gambling_LR_RL/$i/output gambling_LR_RL/$i/output_z tfMRI_GAMBLING $i ev_4dfiles_LR_RL_wrapper
done 


mkdir motor_LR_RL
for i in $motor;
do
	mkdir motor_LR_RL/$i
	python wrapper_LR_RL.py motor_LR_RL/$i/output motor_LR_RL/$i/output_z tfMRI_MOTOR $i ev_4dfiles_LR_RL_wrapper
done

mkdir relational_LR_RL
for i in $relational;
do
	mkdir relational_LR_RL/$i
	python wrapper_LR_RL.py relational_LR_RL/$i/output relational_LR_RL/$i/output_z tfMRI_RELATIONAL $i ev_4dfiles_LR_RL_wrapper
done 
 
