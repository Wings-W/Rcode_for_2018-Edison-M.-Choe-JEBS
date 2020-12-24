I write these codes for my final project.

Recently, I read the article: Optimizing the use of response times for item selection in CAT, which published at JEBS in 2018.
Being interested in the Generalized MIT (GMIT) method, I try to present the results of study 1 with set1 under the "v = 1, w = 0.5" condition.

Three item selection methods were compared (ie., Random selection method, maximum Fisher information selection method and GMIT method with v = 1, w = 0.5).
It should noticed that I didn't implement any repetition for time saving in my presentation.


THE PROJECT FRAMEWORK:

#define function

A. define function1: compute Fisher information for dichotomous IRT model

B. define function2: generate response data for an examinee

C. define function3: generate response time for an examinee

D. define function4: Random selection method

E. define function5: maximum Fisher information selection method

F. define function6: Generalized MIT method

G. define function7: estimate theta for an examinee with maximum likelihood method

H. define function8: estimate theta for an examinee with expected a posteriori method (as an interim substitute of MLE)

I. define function9: estimate speed parameter for an examinee with maximum likelihood method


#set variables

I set some variables used in this study and generate item and person parameters in this part.

#main body

The estimation main body and the computation of seven evaluation criteria (ie., bias and RMSE of theta, bias and RMSE of tau, M and SD of testing times, mean of testing overlap rate).

#save results

Save results and return evaluation indexs and timecost to the console.

You can use these codes to try using response times for item selection in CAT! Also, feel free to ask me any question.

Reference

Choe, E. M., Kern, J. L., & Chang, H. H. (2018). Optimizing the use of response times for item selection in computerized adaptive testing. Journal of Educational and Behavioral Statistics, 43(2), 135-158.
