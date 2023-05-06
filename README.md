# A/B Test for GlowBox
---

## Project Overview

In this project we conduct and analyse an A/B test which will be used for a presentation of data-driven recommendations based on our findings. The A/B test will be conducted in Python.

> An A/B test is an experimentation technique used by businesses to compare two versions of a webpage, advertisement, or product feature to determine which one performs better. By randomly assigning customers or users to either the A or B version, the business can determine which version is more effective at achieving a particular goal.

## Project Background

The e-commerce company Glowbox specializes in sourcing unique and high-quality prodects from around the world. 

The company is primarily known amongst its customer base for boutique fashion items and high-end decor products. However, their food and drink offerings have grown tremendously in the last few months, and the company wants to bring awareness to this product category to increase revenue.

The Growth team decided to run an A/B test that highlights key products in the food and drink category as a banner at the top of the website. 

## Test Groups

**Group A: Control** - existing landing page

**Group B: Treatment** - landing page with food and drink banner

## The setup of the A/B test is as follows:

1. The experiment is only being run on the mobile website.
2. A user visits the GloBox main page and is randomly assigned to either the control or test group. This is the join date for the user.
3. The page loads the banner if the user is assigned to the test group, and does not load the banner if the user is assigned to the control group.
4. The user subsequently may or may not purchase products from the website. It could be on the same day they join the experiment, or days later. If they do make one or more purchases, this is considered a “conversion”.

## The Dataset
GloBox stores its data in a relational database, which you will access through bit.io.

🔗 [Link to database](https://bit.io/griffinmasterschool/ab_test_project)



