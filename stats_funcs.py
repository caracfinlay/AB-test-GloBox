import numpy as np

def get_mean_and_se(df, col):
    '''
    Compute the mean and standard error of a column in a pandas DataFrame.

    Parameters
    ----------
    df : pandas DataFrame
        The DataFrame containing the column to compute the mean and standard error for.
    col : str
        The name of the column to compute the mean and standard error for.

    Returns
    -------
    mean,se: tuple
        A tuple containing two values: the mean of the specified column, and the standard error of the mean.
    '''
    mean = df[col].mean()
    se = df[col].std() / np.sqrt(len(df))
    return mean, se

def critical_value(sig, dof, test_kind='two-tail'):
    '''
    Calculate critical value for various tests.

    Parameters:
    -----------
    sig: float
        Significance level between 0 and 1 e.g. 0.05
    dof: int
        Degrees of freedom i.e. number of data points - 1
    test_kind: string
        Type of test i.e. left, right, two-tail. 
        Default is two-tail. 

    Returns: 
    --------
    critical_value: float 
        The critical value.
    '''
    
    from scipy.stats import t
    if test_kind == 'left':
        crit = t.ppf(q=sig, df=dof)
    elif test_kind == 'right':
        crit = t.ppf(q=1-sig, df=dof)
    else:
        crit = t.ppf(q=1-sig/2, df=dof)
    return crit
def confidence_interval(mean, se, critical_value):
    '''
    Compute a confidence interval for a given mean and standard error.

    Parameters
    ----------
    mean : float
        The sample mean to compute the confidence interval for.
    se : float
        The standard error of the mean.
    critical_value : float
        The critical value for the desired confidence level.

    Returns
    -------
    tuple
        A tuple containing the lower and upper bounds of the confidence interval.
    '''

    lower = mean - critical_value * se
    upper = mean + critical_value * se
    return lower, upper 

def t_test_statistic_1_samp (sample, sample_col, hypothesis_mean):
    '''
    Calculates a t-test for one sample.
       T = (x̄ - μ₀) / (s / √n)
     
    Parameters:
    -----------
    sample: df
        This is you sample dataframe.
    sample_col: string
        This is the column to be aggregated. 
    hypothesis_mean: float
        μ₀: This is the hypothesized population mean, which is the value you want to test against.
    
    Returns:
    --------
    T: float
        This is the calculated t-value, which measures the difference between the sample mean (x̄) 
        and the hypothesized population mean (μ₀) in units of the standard error of the sample mean.
        The resulting t-value is compared to a t-distribution to determine the probability of 
        observing such a large difference between the sample mean and the hypothesized population mean 
        by chance alone. If this probability is low enough (usually set at a significance level of 0.05), 
        we reject the null hypothesis that the sample mean is not significantly different 
        from the hypothesized population mean.
    '''
    sample_n = len(sample)
    mean, se = get_mean_and_se(sample, sample_col)

    t_1_samp = (mean - hypothesis_mean) / se
    return t_1_samp

def t_test_statistic_2_samp (control_sample, treatment_sample, control_col, treatment_col, hypothesis_mean=0):
    '''
    Calculates the t-test statistic for two independent samples.

    Parameters:
    -----------
    sample_1 : pandas DataFrame
        The first DataFrame used for computing the mean and standard error.
    sample_2 : pandas DataFrame
        The second DataFrame used for computing the mean and standard error.
    col_1 : str
        The column used for analysis within the first DataFrame.
    col_2 : str
        The column used for analysis within the second DataFrame. 
    hypothesis_mean : float, optional
        The null hypothesis mean, which is the value being tested against the sample means to 
        determine whether the difference between them is statistically significant. The default is 0. 

    Returns:
    --------
    float
        The t-test statistic for two independent samples.
    '''
    control_mean, control_se = get_mean_and_se(control_sample, control_col)
    treatment_mean, treatment_se = get_mean_and_se(treatment_sample, treatment_col)

    t_2_samp = ((treatment_mean - control_mean) - hypothesis_mean) / np.sqrt((control_se**2 + treatment_se**2))
    return t_2_samp

def confidence_interval_diff_mean(control_sample, treamtent_sample, control_col, treatment_col, sig):
    '''
    Compute a confidence interval for the difference between the means of two samples.

    Parameters
    ----------
    control_sample : pandas DataFrame
        The control sample data containing the column specified by `col_1`.
    treamtent_sample : pandas DataFrame
        The treatment sample data containing the column specified by `col_2`.
    col_1 : str
        The name of the column in `control_sample` to compute the mean and standard error for.
    col_2 : str
        The name of the column in `treatment_sample` to compute the mean and standard error for.
    sig : float
        The desired significance level (e.g., 0.05).

    Returns
    -------
    tuple
        A tuple containing the lower and upper bounds of the confidence interval for the difference between the means of
        the two samples.
    '''
    control_mean, control_se = get_mean_and_se(control_sample, control_col)
    treatment_mean, treamtent_se = get_mean_and_se(treamtent_sample, treatment_col)
    dof = (len(control_sample) - 1) + (len(treamtent_sample) - 1)
    test_kind = 'two-tail'

    sample_stat = treatment_mean - control_mean
    se = np.sqrt(control_se**2 + treamtent_se**2)

    crit_value = critical_value(sig, dof, test_kind='two-tail')

    lower = sample_stat - crit_value * se
    upper = sample_stat + crit_value * se

    return lower, upper

def conversion_rate(df, col):
    '''
    Calculate the conversion rate of a column in a pandas DataFrame.

    Parameters:
    df (pandas.DataFrame): The DataFrame containing the data.
    col (str): The name of the column to calculate the conversion rate for.

    Returns:
    float: The conversion rate as a percentage, i.e., the proportion of rows
           in the DataFrame where the value in the specified column is non-zero.

    Example:
    >>> data = pd.DataFrame({'col1': [0, 1, 0, 1], 'col2': [1, 0, 1, 0]})
    >>> conversion_rate(data, 'col1')
    0.5
    '''

    con_rate = len(df[df[col] != 0]) / len(df)

    return con_rate
def proportion_stats(df, col):

    '''
    Calculate the sample proportion and standard error for a column in a pandas DataFrame.

    Parameters:
    df (pandas.DataFrame): The DataFrame containing the data.
    col (str): The name of the column to calculate the sample proportion and standard error for.

    Returns:
    tuple: A tuple containing the sample proportion and standard error as floats.

    Example:
    >>> data = pd.DataFrame({'col1': [0, 1, 0, 1], 'col2': [1, 0, 1, 0]})
    >>> proportion_stats(data, 'col1')
    (0.5, 0.2886751345948129)
    '''

    proportion = len(df[df[col] != 0]) / len(df)
    se = np.sqrt(proportion*(1-proportion)/len(df))
    
    return proportion, se
def critical_value_proportion(sig, test_kind='two-tail'):

    '''Calculate critical value for various tests.

    Parameters:
    -----------
    sig: float
        Significance level between 0 and 1 e.g. 0.05
    test_kind: string
        Type of test i.e. left, right, two-tail. 
        Default is two-tail. 

    Returns: 
    --------
    critical_value: float 
        The critical value for the given significance level and test type.
        
    Example:
    >>> critical_value_proportion(0.05, 'two-tail')
    1.959963984540054
        
    '''
    from scipy.stats import norm
    if test_kind == 'left':
        crit = norm.ppf(q=sig)
    elif test_kind == 'right':
        crit = norm.ppf(q=1-sig)
    else:
        crit = norm.ppf(q=1-sig/2)
    return crit
def pooled_prop_se(p1, p2, df1, df2):
    
    n1 = len(df1)
    n2 = len(df2)

    p = (p1 * n1 + p2 * n2) / (n1 + n2)

    pooled_se =  np.sqrt(p*(1-p) * (1/n1 + 1/n2))

    return pooled_se
def pooled_confidence(sample_stat, critical_value_proportion, pooled_prop_se):

    lower = sample_stat - critical_value_proportion * pooled_prop_se
    upper = sample_stat + critical_value_proportion * pooled_prop_se

    return lower, upper

def pooled_t_statistic_prop(p1 , p2, df1, df2 ):

    '''Computes the pooled t-statistic for two independent samples with known proportions.

    Parameters
    ----------
    p1 : float
        Proportion of successes in sample 1.
    p2 : float
        Proportion of successes in sample 2.
    df1 : array-like
        Sample 1.
    df2 : array-like
        Sample 2.

    Returns
    -------
    pooled_t : float
        The pooled t-statistic.

    Notes
    -----
    The pooled t-statistic is used to test the hypothesis that two independent samples have equal population means. It assumes that the variances of the two populations are equal. The formula for the pooled t-statistic is:
    pooled_t = (p1 - p2) / sqrt(p * (1 - p) * (1 / n1 + 1 / n2))

        where:
            - p is the pooled proportion of successes
            - n1 and n2 are the sample sizes of sample 1 and sample 2, respectively'''


    n1 = len(df1)
    n2 = len(df2)

    p = (p1 * n1 + p2 * n2) / (n1 + n2)

    pooled_t = (p1-p2) / np.sqrt( p*(1-p) * (1/n1 + 1/n2) )

    return pooled_t
def se_unpooled_prop(p1 , p2, df1, df2):
    '''Computes the unpooled standard-error for z-interval for two independent samples with known proportions.

    Parameters
    ----------
    p1 : float
        Proportion of successes in sample 1.
    p2 : float
        Proportion of successes in sample 2.
    df1 : array-like
        Sample 1.
    df2 : array-like
        Sample 2.

    Returns
    -------
    unpooled_se : float
        The unpooled standard-error.'''
    
    n1 = len(df1)
    n2 = len(df2)
    
    unpooled_se = np.sqrt((p1*(1-p1)/n1) + (p2*(1-p2)/n2))

    return unpooled_se
# incomplete

def unpooled_t_statistic_prop(p1 , p2, df1, df2 ):

    '''Computes the unpooled t-statistic for two independent samples with known proportions.

    Parameters
    ----------
    p1 : float
        Proportion of successes in sample 1.
    p2 : float
        Proportion of successes in sample 2.
    df1 : array-like
        Sample 1.
    df2 : array-like
        Sample 2.

    Returns
    -------
    pooled_t : float
        The pooled t-statistic.'''


    n1 = len(df1)
    n2 = len(df2)

    

    return unpooled_t

def determine_pool( h0):
    if h0 == 0:
        print("The proportions are equal, therefore use the pooled standard error. ")
    else:
        print("The proportions are not equal, therefore use the unpooled standard error.")