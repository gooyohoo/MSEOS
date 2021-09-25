# -*- coding: utf-8 -*-
"""
Created on Sun Mar 12 21:16:05 2017

@author: Huhu
"""
"""""""""""""""""""dateset input"""""""""""""""""""""""""""""""
from os import listdir
files = listdir('.')
if 'blackouts.txt' not in files:
    import urllib as url
    url.request.urlretrieve('https://raw.github.com/jeffalstott/powerlaw/master/manuscript/blackouts.txt', 'blackouts.txt')
if 'words.txt' not in files:
    import urllib as url
    url.request.urlretrieve('https://raw.github.com/jeffalstott/powerlaw/master/manuscript/words.txt', 'words.txt')
if 'worm.txt' not in files:
    import urllib as url
    url.request.urlretrieve('https://raw.github.com/jeffalstott/powerlaw/master/manuscript/worm.txt', 'worm.txt')

from numpy import genfromtxt
blackouts = genfromtxt('blackouts.txt')#/10**3
words = genfromtxt('words.txt')
worm = genfromtxt('worm.txt')
worm = worm[worm>0]
critical_branch=genfromtxt('critical_branch.txt')

"""""""""""""""""""""""figure 1 prepare"""""""""""""""""""""""""""  
import matplotlib.pyplot as plt
import powerlaw as pwl

import pylab
pylab.rcParams['xtick.major.pad']='8'
pylab.rcParams['ytick.major.pad']='8'
#pylab.rcParams['font.sans-serif']='Arial'

from matplotlib import rc
rc('font', family='sans-serif')
rc('font', size=10.0)
rc('text', usetex=False)


from matplotlib.font_manager import FontProperties

panel_label_font = FontProperties().copy()
panel_label_font.set_weight("bold")
panel_label_font.set_size(12.0)
panel_label_font.set_family("sans-serif")

def plot_basics(data, data_inst, fig, units):
    from powerlaw import plot_pdf, Fit, pdf
    annotate_coord = (-.4, .95)
    ax1 = fig.add_subplot(n_graphs,n_data,data_inst)
    x, y = pdf(data, linear_bins=True)#线性的bins,or 对数的bins
    ind = y>0#index of y>0
    y = y[ind]
    x = x[:-1]
    x = x[ind]
    ax1.scatter(x, y, color='r', s=.5)
    plot_pdf(data[data>0], ax=ax1, color='b', linewidth=2)
    from pylab import setp
    setp( ax1.get_xticklabels(), visible=False)

    if data_inst==1:
        ax1.annotate("A", annotate_coord, xycoords="axes fraction", fontproperties=panel_label_font)
    
    from mpl_toolkits.axes_grid.inset_locator import inset_axes
    ax1in = inset_axes(ax1, width = "30%", height = "30%", loc=3)
    ax1in.hist(data, normed=True, color='b')
    ax1in.set_xticks([])
    ax1in.set_yticks([])

    
    ax2 = fig.add_subplot(n_graphs,n_data,n_data+data_inst, sharex=ax1)
    plot_pdf(data, ax=ax2, color='b', linewidth=2)
    fit = Fit(data, xmin=1, discrete=True)#设置了最小值
    fit.power_law.plot_pdf(ax=ax2, linestyle=':', color='g')
    #p = fit.power_law.pdf()

    ax2.set_xlim(ax1.get_xlim())
    
    fit = Fit(data, discrete=True)
    fit.power_law.plot_pdf(ax=ax2, linestyle='--', color='g')
    from pylab import setp
    setp( ax2.get_xticklabels(), visible=False)

    if data_inst==1:
       ax2.annotate("B", annotate_coord, xycoords="axes fraction", fontproperties=panel_label_font)        
       ax2.set_ylabel(u"p(X)")# (10^n)")
        
    ax3 = fig.add_subplot(n_graphs,n_data,n_data*2+data_inst)#, sharex=ax1)#, sharey=ax2)
    fit.power_law.plot_pdf(ax=ax3, linestyle='--', color='g')
    fit.exponential.plot_pdf(ax=ax3, linestyle='--', color='r')#指数分布的概率密度函数
    fit.plot_pdf(ax=ax3, color='b', linewidth=2)
    
    ax3.set_ylim(ax2.get_ylim())
    ax3.set_xlim(ax1.get_xlim())
    
    if data_inst==1:
        ax3.annotate("C", annotate_coord, xycoords="axes fraction", fontproperties=panel_label_font)

    ax3.set_xlabel(units)
    
"""""""""""""""""""""""figure 1"""""""""""""""""""""""""""    
n_data = 3
n_graphs = 4
f = plt.figure(figsize=(8,11))

data = words
data_inst = 1
units = 'Word Frequency'
plot_basics(data, data_inst, f, units)

data_inst = 2
data =critical_branch #worm
units = 'Neuron Connections'
plot_basics(data, data_inst, f, units)

data = blackouts
data_inst = 3
units = 'Population Affected\nby Blackouts'
plot_basics(data, data_inst, f, units)

f.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=.3, hspace=.2)
figname = 'FigWorkflow'
f.savefig(figname+'.eps', bbox_inches='tight')
#f.savefig(figname+'.tiff', bbox_inches='tight', dpi=300)

"""""""""""""""""""""detail"""""""""""""""""""""""""""""
"http://nbviewer.jupyter.org/github/jeffalstott/powerlaw/blob/master/manuscript/Manuscript_Code.ipynb"
#数据拟合和数据分布匹配
data=words
fit=pwl.Fit(data,discrete=True)
fit.alpha
fit.sigma#std
fit.distribution_compare('power_law', 'exponential')#数据分布匹配，返回配对结果和p值

plt.figure()
pwl.plot_pdf(data,color='b')#对数bins的拟合
pwl.plot_pdf(data, linear_bins = True, color = 'r')
fit.distribution_compare('power_law', 'exponential')
#数据可视化，pdf,cdf,ccdf
plt.figure()
data=words
fit=pwl.Fit(data,discrete=True)
fig2=fit.plot_pdf(color='r',linewidth=2)#数据平滑后的pdf
fit.power_law.plot_pdf(ax=fig2,color = 'r', linestyle = '--')#拟合的最接近的powerlaw pdf
fit.plot_ccdf(ax=fig2,color = 'b')#补cdf
fit.power_law.plot_ccdf(ax=fig2,color = 'b', linestyle = '--')#拟合的powerlaw pdf

x,y=fit.cdf()#fit 后的累积概率密度函数，参数选择可以设置为fit前的data
bin_edges,probability=fit.pdf()
fit.lognormal.cdf(data=[300,350])
y=fit.lognormal.pdf()#对数正太分布的概率密度函数

#指定最小fit点和搜寻最小fit点 Identifying the Scaling Range
fit=pwl.Fit(blackouts)
fit.xmin
fit.fixed_xmin
fit.power_law.alpha#结果是一样的
fit.alpha#结果是一样的
fit.power_law.D#拟合和真实分布间的距离Kolmogorov-Smirnov distance
#####################################################
fit=pwl.Fit(blackouts,xmin=1.0)#指定最小拟合开始点
fit.xmin
fit.fixed_xmin
fit.power_law.alpha#结果是一样的
fit.alpha#结果是一样的
fit.power_law.D#拟合和真实分布间的距离Kolmogorov-Smirnov distance

fit=pwl.Fit(words,xmin = (250.0, 300.0))
fit.xmin
fit.fixed_xmin
fit.given_xmin
#指定最大fit点 Identifying the Scaling Range
fit=pwl.Fit(words,xmax=1000.0)

"""""""""""""""""""""""figure 3"""""""""""""""""""""""""""
data = words
plt.figure()##ccdf容易看出xmin. 指数分布的cdf任然为等指数的幂律分布
fit = pwl.Fit(data, discrete=True, xmax=None)
FigCCDFmax = fit.plot_ccdf(color='b', label=r"Empirical, no $x_{max}$")
fit.power_law.plot_ccdf(color='b', linestyle='--', ax=FigCCDFmax, label=r"Fit, no $x_{max}$")
fit = pwl.Fit(data, discrete=True, xmax=1000)
fit.plot_ccdf(color='r', label=r"Empirical, $x_{max}=1000$")
fit.power_law.plot_ccdf(color='r', linestyle='--', ax=FigCCDFmax, label=r"Fit, $x_{max}=1000$")

FigCCDFmax.set_ylabel(u"p(X≥x)")
FigCCDFmax.set_xlabel(r"Word Frequency")
handles, labels = FigCCDFmax.get_legend_handles_labels()
leg = FigCCDFmax.legend(handles, labels, loc=3)
leg.draw_frame(False)

figname = 'FigCCDFmax'
plt.savefig(figname+'.eps', bbox_inches='tight')

#分布的参数
fit.power_law
fit.power_law.parameter1_name
fit.power_law.parameter1

fit.lognormal.mu
fit.lognormal.parameter1_name
fit.lognormal.parameter2_name
fit.lognormal.parameter3_name==None
#分布比较，确定powerlaw is the best description of it
R, p = fit.distribution_compare('power_law', 'exponential',normalized_ratio = True)
print( R,p)
fit.distribution_compare('power_law', 'truncated_power_law')#truncated_power_law is better with a size bound 

"""""""""""""""""""""""figure 4"""""""""""""""""""""""""""
data = words
plt.figure()
fit = pwl.Fit(data, discrete=True)
fit.distribution_compare('power_law','lognormal')
fig4=fit.plot_ccdf(linewidth=3)
fit.power_law.plot_ccdf(ax=fig4,color='r',linestyle='--')
fit.lognormal.plot_ccdf(ax=fig4,color='g',linestyle='--')

#create simulated data to validate fit quality
data = blackouts
fit = pwl.Fit(data)
simulated_data = fit.power_law.generate_random(10000)#way1


theoretical_distribution = pwl.Power_Law(xmin =5.0,parameters =[2.5])#way2
simulated_data = theoretical_distribution.generate_random(10000)
fit=pwl.Fit(simulated_data)
fit.power_law.xmin,fit.power_law.alpha

plt.figure()
pwl.plot_pdf(simulated_data,linewidth=3)
fit.power_law.plot_pdf(simulated_data,linestyle='--',color='r')

#restricted parameter range
data=words
fit=pwl.Fit(data,sigma_threshold=.1)
fit.alpha,fit.sigma,fit.xmin

parameter_range={'alpha':[2.3,None],'sigma':[None,.2]}#None is necessary
fit=pwl.Fit(data,parameter_range=parameter_range)
fit.alpha,fit.sigma,fit.xmin

parameter_range = lambda self: self.sigma/self.alpha < .5 #None is necessary
fit=pwl.Fit(data,parameter_range=parameter_range)
fit.alpha,fit.sigma,fit.xmin

#Multiple Possible Fits choose
data = blackouts
fit = pwl.Fit(data, sigma_threshold=.1)
print(fit.xmin, fit.D, fit.alpha)

fit = pwl.Fit(data)
print(fit.xmin, fit.D, fit.alpha)
plt.figure()
from matplotlib.pylab import plot
plot(fit.xmins, fit.Ds, label=r'$D$')
plot(fit.xmins, fit.sigmas, label=r'$\sigma$', linestyle='--')
plot(fit.xmins, fit.sigmas/fit.alphas, label=r'$\sigma /\alpha$', linestyle='--')
####
plt.ylim(0, .4)
plt.legend(loc=4)
plt.xlabel(r'$x_{min}$')
plt.ylabel(r'$D,\sigma,\alpha$')

#Selecting xmin with Other Distance Metrics
fit = pwl.Fit(data, xmin_distance = 'D')
fit = pwl.Fit(data, xmin_distance = 'V')
fit = pwl.Fit(data, xmin_distance = 'Asquare')





