#defines the class of the polar plot
class polar_axis():
  step=5
  def __init__(self,fig,ax,data):
    self.fig = fig #pylab.figure()
    self.ax = ax #pylab.axes([0.025,0.025,0.95,0.95], polar=True)
    self.data=data
    self.reset=True
    self.ax.cla()
    datatmp=self.data
    self.subdata=self.dump_timing_level(datatmp)
    self.draw_polarplot()
    self.fig.canvas.mpl_connect('motion_notify_event',self._info_callback)
    #self.fig.canvas.mpl_connect('pick_event',self._new_level)
    self.fig.canvas.mpl_connect('button_press_event',self._new_level)

  def draw_polarplot(self):
    import pylab
    import numpy as np
    #self.subdata=data
    self.tot=self.subdata["time"][0]
    self.times=self.subdata["time"]
    self.N=len(self.times)

    self.width=np.array(self.subdata["time"])/self.tot*2*np.pi
    self.theta=np.array(self.subdata["theta"])/self.tot*2*np.pi
    self.bot=np.array(self.subdata["level"])
    self.names=self.subdata["names"]
    self.radii = np.array(self.N*[self.step])

    self.bars = self.ax.bar(self.theta,self.radii,
                          width=self.width,
                          bottom=self.step*self.bot,picker=True)

    ilev=0
    #maxlev=max(self.bot)
    for r,bar,ilev in zip(self.radii, self.bars,self.theta):
       #print ilev,'hello',float(ilev)/float(N),maxlev
       #bar.set_facecolor( pylab.cm.jet(float(ilev)/maxlev))
       bar.set_facecolor(pylab.cm.jet(float(ilev)/(2.0*np.pi)))
       bar.set_alpha(0.5)
       ilev+=1

    self.ax.set_xticklabels([])
    self.ax.set_yticklabels([])
    self._wrotesome=False
    self.info=self.ax.text(0,0,'',fontsize=150)

  ## to be used what it is the moment to show the plot(s)
  #def show(self):
  #  import pylab
  #  try:
  #    pylab.show()
  #  except KeyboardInterrupt:
  #    raise

  def _find_name(self,th,level):
    import numpy as np
    levth=[]
    ipiv=[]
    for i in range(self.N):
      if self.bot[i]==level:
        levth.append(self.theta[i])
        ipiv.append(i)
    if len(levth)==0: return None
    to_find_zero=np.array(levth)-th
    to_find_zero[to_find_zero > 0]=-4*np.pi
    routine=np.argmin(-to_find_zero)
    if routine == len(to_find_zero)-1 and abs(routine)> 0.5:
      return None
    else:
      return (self.theta[ipiv[routine]],self.names[level][routine],self.times[ipiv[routine]])
    #if min(abs(to_find_zero))==0.0:
    #  routine=np.argmin(abs(to_find_zero))
    #  return (th,self.names[level][routine],self.times[ipiv[routine]])
    #else:
    #  to_find_zero[to_find_zero > 0]=-4*np.pi
    #  routine=np.argmin(-to_find_zero)
    #  return (self.theta[ipiv[routine]],self.names[level][routine],self.times[ipiv[routine]])


  def _info_string(self,xdata,level):
    #create string to be plotted in plot
    (tt,name,time)=self._find_name(xdata,level)
    info=name+':\n'+str(time)+" s ("
    info+= " %6.2f %s" % (time/self.tot*100.0,"% ) \n")
    #find lower level
    it=range(level)
    it.reverse()
    parents=''
    for i in it:
      (tt,name,time)=self._find_name(tt,i)
      parents+="-"+name+"("+str(time)+"s,"
      parents+= " %6.2f %s" % (time/self.tot*100.0,"% ) \n")
    if len(parents) > 0:
      info+="Calling path:\n"+parents
    return info

  def _new_level(self,event):
    #click=event.mouseevent.button
    click=event.button
    datatmp=self.data
    if click==3:
      if not self.reset:
        self.ax.texts.remove(self.comebackinfo)
        self.reset=True
      nameandlev=None
    elif click==1:
      #thisline = event.artist
      #xdata, ydata = thisline.get_xy()
      xdata, ydata=event.xdata,event.ydata
      level = int(ydata)/self.step
      goodpoint=self._find_name(xdata,level)
      if goodpoint is None: return
      (tt,name,time)=goodpoint
      if not self.reset: return #no double level
      self.reset=False
      nameandlev=(name,level)
      #extract the data to be plotted for the internal level
    else:
      return
    self.subdata=self.dump_timing_level(datatmp,starting_point=nameandlev)
    self.ax.cla()
    self.draw_polarplot()
    if not self.reset:
      self.comebackinfo=self.ax.text(0.02,0.95,
                                     'Subroutines of '+name+' (Right click to reset)',
                                     fontsize=15,transform = self.ax.transAxes)
    self.fig.canvas.draw()

  def _info_callback(self,event):
    #print dir(event)
    #thisline = event.artist
    #xdata, ydata = thisline.get_xy()
    xdata=event.xdata#data
    ydata=event.ydata#data
    #print 'test',xdata,ydata,event.x,event.y
    if xdata is None or ydata is None: return
    level = int(ydata)/self.step
    goodpoint=self._find_name(xdata,level)

    #then plot the string
    if self._wrotesome:
      self.ax.texts.remove(self.info)
      self.info=self.ax.text(0,0,'',fontsize=150)
      self._wrotesome=False
    if goodpoint is None: return
    offset = 0.02
    #print 'nowinfo',info
    self.info= self.ax.text( offset, offset, self._info_string(xdata,level),
                             fontsize = 15,transform = self.ax.transAxes )
    self._wrotesome=True
    self.fig.canvas.draw()

  def dump_timing_level(self,level,starting_point=None,ilev=0,theta=0,data=None): ##data={"time":[],"level":[],"theta":[],"names":[[]]},reset=False):
    """Inspect the first level of the given dictionary and dump the profile subroutines at this level"""
    import pylab
    if data is None: data={"time":[],"level":[],"theta":[],"names":[[]]}
    subs=data["time"]
    tht=data["theta"]
    lev=data["level"]
    nms=data["names"]
    if ilev == len(nms):
      nms.append([])
    #tel=0
    tel=theta #entry point of the routine
    t0=0
    start=starting_point
    import copy
    leveltmp=copy.deepcopy(level)
    brk=False
    for routine in leveltmp:
      #first eliminate the subroutines from the level
      try:
        sublevel=routine.pop("Subroutines")
      except:
        sublevel=None
      for name in routine:
        if starting_point is not None:
          if (name,ilev) != starting_point:
            continue
          else:
            brk=True
            start=None #from now on take the data
            ilev=0
        t0=routine[name][0] #take the time
        subs.append(t0) #take the time
        tht.append(tel)
        lev.append(ilev)
        nms[ilev].append(name)
      if sublevel is not None:
        jlev=ilev+1
        self.dump_timing_level(sublevel,starting_point=start,ilev=jlev,theta=tel,data=data)
      tel=tel+t0
      if brk: break
    return data

class TimeData:
  barwidth=0.9#0.35
  ignored_counters=['CPU parallelism','Routines timing and number of calls','SUMMARY','Report timestamp','Hostnames']
  def __init__(self,*filenames,**kwargs):
    """
    Class to analyze the timing from a Futile run
    ------
    Arguments
    filenames: list of the files which have to be treated
    plottype: Decide the default yscale for the plotting
              May be "Percent" or "Seconds"
    static: Show the plot statically for screenshot use
    fontsize: Determine fontsize of the bar chart plot
    nokey: Remove the visualization of the key from the main plot
    counter: retrieve as reference value the counter provided
    only_last: retrieve only the last document in the time.yaml
    """
    #here a try-catch section should be added for multiple documents
    #if (len(filename) > 1
    from futile import Figures as Fig
    only_last = kwargs.get('only_last',False)
    self.log=[]
    for filename in filenames:
      from futile import YamlIO as Yaml
      tmplogs=Yaml.load(filename,doc_lists=True)
      if only_last:
        self.log+=[[ a for a in tmplogs][-1]]
      else:
        self.log+=tmplogs
      #try:
      #  self.log+=[yaml.load(open(filename, "r").read(), Loader = yaml.CLoader)]
      #except:
      #  self.log+=yaml.load_all(open(filename, "r").read(), Loader = yaml.CLoader)
    #create the figure environemnt
    self.figures=Fig.FigureSet(title='Profiling',QuitButton=True,twinaxes=True)
    #self.barfig = None
    #self.axbars = None
    #self.newfigs =[]
    self.radio = None
    self.toggle_unbalancing = False
    self.quitButton = None
    self.lined=[]
    self.plot_start=kwargs.get('plottype','Seconds')
    self.static = kwargs.get('static',False)
    self.fontsize=kwargs.get('fontsize',15)
    self.nokey=kwargs.get('nokey',False)
    self.strong_scaling=kwargs.get('strong_scaling',True)
    counter=kwargs.get('counter','WFN_OPT') #the default value, to be customized
    if counter in self.counters():
      self.inspect_counter(counter)
    else:
      print("Warning: counter not initialized, check the available counters")

  def counters(self):
    "Inspect the available counters"
    cnts=[]
    for doc in self.log:
      if doc is None: continue
      for internal_cnts in doc.keys():
        if internal_cnts not in cnts and internal_cnts not in self.ignored_counters: cnts.append(internal_cnts)
    return cnts

  def _refill_classes(self,main):
    self.classes=[]
    self.totals=[]
    for doc in self.log:
        if doc is None: continue
        scf=doc.get(main)
        if scf is None: continue
        self.scf.append(scf)
        self.routines.append(doc.get("Routines timing and number of calls"))
        self.hostnames.append(doc.get("Hostnames"))
        self.unbalancings.append(True)
        loclass=scf["Classes"].keys()
        if 'Total' in loclass:
          self.totals.append(scf["Classes"]["Total"][1])
        for cs in loclass:
          if cs not in self.classes and cs != "Total": self.classes.append(cs)
          if len(scf["Classes"][cs]) == 2: self.unbalancings[-1]=False
        mpit=doc.get("CPU parallelism")
        title=doc.get("Run name","Unknown")
        if mpit is not None:
            mpi=mpit.get('MPI tasks')
            omp=mpit.get('OMP threads')
            if title is "Unknown":
              title=str(mpi) if omp is None else str(mpi)+'-'+str(omp)
            ncores=mpi
            if omp is not None: ncores*=omp
        else:
            ncores=0.0
        self.ids.append(title)
        self.ncores.append(ncores)
    self.classes.sort()
    self.classes.append("Unknown") #an evergreen

  def inspect_counter(self,counter,unit=None):
    self.routines=[]
    self.hostnames=[]
    self.unbalancings=[]
    self.scf=[]
    self.ids=[]
    self.ncores=[] #number of cores for each of the valid run
    self.vals=self.plot_start if unit is None else unit
    self._refill_classes(counter)
    #here we might add the policy to add new figures or delete previous ones
    data=self.collect_categories(self.scf,self.vals)
    fig,axis=self.figures[0]
    axis.cla()
    self.draw_barfigure(fig,axis,data,title="Counter "+counter)

  def inspect_category(self,cat):
    data=self.find_items(cat,self.scf)
    if self.figures.exists(cat):
      newfig,newax=self.figures[cat]
      #erase the plot in case it is already available
      newax.cla()
    else:
      newfig,newax=self.figures.add(title=cat,twinaxes=True)
    self.draw_barfigure(newfig,newax,data,title=cat)
    newfig.show()

  def show(self):
    self.figures.show()

  def _strong_scaling_lines(self,data,axis):
    totals=data[1]
    ref=totals[0] # the reference value is assumed to be the first
    coreref=float(self.ncores[0])
    #print self.totals
    if self.vals=='Seconds':
      label='Speedup'
      dataline=[ref/b for b in totals]
    else:
      label='Parallel Efficiency'
      dataline=[ (ref/b)/(c/coreref)  for b,c in zip(totals,self.ncores)]
    self.draw_lineplot(axis.twin,dataline,label)


  def draw_barfigure(self,fig,axis,data,title):
    import numpy as np
    from matplotlib.widgets import RadioButtons,Button
    import matplotlib.pyplot as plt
    from futile import Figures as Fig
    if self.static: fig.patch.set_facecolor("white")
    #write the plot in the axis
    self._draw_barplot(axis,data[0],self.vals,title=title,nokey=self.nokey)
    #add the lines useful for a strong scaling test
    if self.strong_scaling: self._strong_scaling_lines(data,axis)
    #if self.vals == 'Percent': axis.set_yticks(np.arange(0,100,10))
    #then add buttons to the plot
    if self.radio is None and not self.static:
      try:
          self.radio = RadioButtons(plt.axes([0.0, 0.75, 0.08, 0.11], axisbg='lightgoldenrodyellow'), ('Seconds', 'Percent'),
                                    active=1 if self.vals=='Percent' else 0)
      except AttributeError:
          self.radio = RadioButtons(plt.axes([0.0, 0.75, 0.08, 0.11], facecolor='lightgoldenrodyellow'), ('Seconds', 'Percent'),
                                    active=1 if self.vals=='Percent' else 0)
      self.radio.on_clicked(self.replot)
      tt=axis.get_xticks()
      routics=[Fig.axis_from_data(fig,axis,(tic-0.45*self.barwidth,self.barwidth)) for tic in tt]
      unbtics=[Fig.axis_from_data(fig,axis,(tic+0.15*self.barwidth,self.barwidth)) for tic in tt]
      self.routine_buttons=[]
      from functools import partial
      for i,t in enumerate(routics):
        if self.routines[i] is None: continue
        but=Button(plt.axes([t[0], 0.0, 0.15*t[1], 0.05]), 'R')
        but.on_clicked(partial(self.routines_plot,i))
        self.routine_buttons.append(but)
      for i,t in enumerate(unbtics):
        if not self.unbalancings[i]: continue
        but=Button(plt.axes([t[0], 0.0, 0.15*t[1], 0.05]), 'U')
        but.on_clicked(partial(self.workload_plot,i))
        self.routine_buttons.append(but)
      #fig.canvas.mpl_connect('pick_event',self._onclick_ev)
      fig.canvas.mpl_connect('button_press_event',
                             partial(self._onclick_ev,fig,axis))
    fig.canvas.draw()

  def routines_plot(self,index,event=None):
    "Draw the plot of the routines level for the run identified by index"
    toplt=self.routines[index]
    title='Routines '+str(self.ids[index])+'-'+str(index)
    if self.figures.exists(title): return
    #data=dump_timing_level(toplt)#,starting_point='cluster')
    fig,ax=self.figures.add(title=title,
                            axes={'rect':[0.025,0.025,0.95,0.95],
                                  'polar':True})
    if not hasattr(self,'routine_plots'): self.routine_plots=[]
    plt=polar_axis(fig,ax,toplt)
    self.routine_plots.append(plt)
    fig.show()

  def workload_plot(self,index,event=None):
    "Draw the plot of the workload of different classes for the run identified by index"
    from functools import partial
    title='Unbalancing '+str(self.ids[index])+'-'+str(index)
    if self.figures.exists(title): return
    fig,ax=self.figures.add(title=title)
    leg,plts,show=self.load_unbalancing(ax,self.scf[index]["Classes"],self.hostnames[index])
    idx=len(self.lined)
    self.lined.append({})
    for legline, origline,yes in zip(leg.get_patches(), plts, show):
          legline.set_picker(5)  # 5 pts tolerance
          self.lined[idx][legline] = origline
          if not yes: self._toggle_visible(origline,legline)
    fig.canvas.mpl_connect('pick_event', partial(self._onpick_curve,idx,fig))
    fig.show()

  def find_items(self,category,dict_list):
    """For a given category find the items which have them"""
    import numpy as np
    items={}
    for idoc in range( len(dict_list) ):
        for cat in dict_list[idoc]["Categories"]:
            dicat=dict_list[idoc]["Categories"][cat]
            if dicat["Class"] == category:
                if cat not in items:
                    items[cat]=np.zeros(len(dict_list))
                items[cat][idoc]=dicat["Data"][self.iprc]
    return [ (cat,items[cat]) for cat in items],[ doc["Classes"][category][1] for doc in dict_list]

  def collect_categories(self,dict_list,vals):
    """ Collect all the categories which belong to all the dictionaries """
    import numpy as np
    if vals == 'Percent':
      self.iprc=0
    elif vals == 'Seconds':
      self.iprc=1
    catsdats=[]
    self.values_legend=[]
    totals=[]
    for cat in self.classes:
      try:
        if cat == "Unknown": #this is always done
            data_unk=[]
            for doc in dict_list:
                totals.append(doc["Classes"]["Total"][1])
                percent_unk=100.0-doc["Classes"]["Total"][0]
                if self.iprc==0:
                    data_unk.append(percent_unk)
                elif self.iprc==1:
                    time_unk=(doc["Classes"]["Total"][1])*percent_unk/100.0
                    data_unk.append(time_unk)
            dat=np.array(data_unk)
        else:
            dat=np.array([doc["Classes"][cat][self.iprc] for doc in dict_list])
        catsdats.append((cat,dat))
        self.values_legend.append(cat)
      except Exception as e:
        print('EXCEPTION FOUND',e)
        print("category",cat,"not present everywhere")
    return catsdats,totals

  def _onclick_ev(self,fig,axis,event):
    import matplotlib.pyplot as plt
    #thisline = event.artist
    #xdata, ydata = thisline.get_xy()
    #print dir(event)
    #print event.inaxes
    xdata,ydata=event.x,event.y #event.xdata,event.ydata
    if xdata is None or ydata is None: return
    #print 'picked',xdata,ydata
    xdata,ydata=axis.transData.inverted().transform((xdata,ydata))
    #print 'transformed',xdata,ydata
    #data_from_data(fig,dst=axis,src=axis.twin,data=(int(xdata),ydata))
    #print 'converted',xdata,ydata
    if ydata < 0.0 or xdata < 0.0 or xdata >= len(self.scf): return
    xdata=int(xdata)
    #find the category which has been identified
    y0data=0.0
    category=None
    for cat in self.values_legend:
      if cat == 'Unknown': continue
      y0data+=self.scf[xdata]["Classes"][cat][self.iprc]
      #print 'cat,y0data',cat,y0data,ydata
      if y0data > ydata:
        category=cat
        break
    #print 'category',category
    if category is not None and category != "Unknown": self.inspect_category(category)
  def _barplot_gnuplot(self,data,ids,lookup=None,select_category=None):
      import Gnuplot
      import numpy as np
      gpdata=[]
      for cat,dat in data:
        if select_category is not None and cat not in select_category: continue
        #gpdata.append(Gnuplot.Data(np.arange(len(dat)),dat,title=cat))
        lup=np.arange(len(dat)) if lookup is None else lookup
        toplt=[dat[0]]+dat[lup].tolist()
        gpdata.append(Gnuplot.Data(toplt,title=cat))
      names=[ids[i] for i in lup]
      return gpdata,names

  def _aggregate_names(self,data,list_agg):
    """
    Aggregate the names of the plot in some different categories
    the structure of dict_agg should ba a list of tuples
    (newkey ,[list of oldkeys])
    between the original keys and the desired keys
    """
    newdata={newkey: 0.0 for newkey,oldkeys in list_agg}
    for cat,dat in data:
      found=False
      for newkey,oldkeys in list_agg:
        if cat in oldkeys:
          newdata[newkey]+=dat
          found=True
          break
      if not found: newdata[cat]=dat
    result=[]
    for newkeys,oldkeys in list_agg:
      result.append((newkeys,newdata.pop(newkeys)))
    for newkeys in newdata:
      result.append((newkeys,newdata[newkeys]))
    return result

  def gnuplot_figure(self,lookup=None,epsfile=None,aggregate=None,select_category=None):
      """Create a figure to be plotted with gnuplot

      Create a gnuplot histogram that can be plotted for production results.

      Args:
         lookup (list): A list of the items of the `TimeData` instance that will be considered for the histogram
         epsfile (str): Name of the eps file in which the data will be plot
         aggregate (list): list of tuples of the categories to aggregate, accoding to the convention
              (newkey ,[list of oldkeys]) between the original keys and the desired keys
         select_category (str): Name of the category to be plot, if specified the internal items of the category will be plot.
      """
      import Gnuplot,numpy as np
      data=self.collect_categories(self.scf,self.vals)
      dataplot=data[0] if aggregate is None else self._aggregate_names(data[0],aggregate)
      gpdata,names=self._barplot_gnuplot(dataplot,self.ids,lookup,select_category)
      gp=Gnuplot.Gnuplot()
      gp("""
set style fill  pattern 1 border -1
set style rectangle back fc  lt -3 fillstyle  solid 1.00 border -1
set key title ""
set key outside right bottom vertical Left reverse enhanced autotitles columnhead nobox
set key invert samplen 4 spacing 1 width 0 height 0
set style increment default
unset style line
unset style arrow
set style histogram rowstacked title  offset character 0, 0, 0
set pointsize 1
set style data histograms
set style function lines
set ytics border in scale 1,0.5 nomirror norotate  offset character 0, 0, 0
set ytics autofreq
set y2tics border in scale 1,0.5 nomirror norotate  offset character 0, 0, 0
set y2tics autofreq
set xtics 0,0.5,0
""")
      gp('set xrange [-1:'+str(len(names))+']')
      if self.vals=='Percent': gp('set yrange [0:100]')
      #gp('set xtics add (" " -1)')
      for i,xtic in enumerate(names):
          gp('set xtics add ("'+str(xtic)+'" '+str(i)+') rotate by 45 right')
      #gp('set xtics add (" " '+str(len(names))+')')
      gp.ylabel(self.vals)
      gp.plot(*gpdata)
      if epsfile is not None: gp.hardcopy(filename=epsfile,eps=True,enhanced=True,color=True)
      return gp

  def _draw_barplot(self,axbars,data,vals,title='Time bar chart',static=False,nokey=False):
    import numpy as np
    from pylab import cm as cm
    ndata=len(data[0][1])
    ind=np.arange(ndata)
    bot=np.array(ndata*[0.0])
    width=self.barwidth
    icol=1.0
    for cat,dat in data:
      #print 'cat',cat,dat
      #for i in range(len(self.ids)):
      #  print self.ids[i],dat[i]
      plt=axbars.bar(ind,dat,width,bottom=bot,color=cm.jet(icol/len(data)),picker=True,label=cat)
      #self.plts.append(plt)
      bot+=dat
      icol+=1.0
    #drawn_classes=np.array(self.values_legend)
    axbars.set_title(title,fontsize=self.fontsize)
    axbars.set_ylabel(vals,fontsize=self.fontsize)
    axbars.set_xticks(ind+width/2.)
    axbars.set_xticklabels(np.array(self.ids),size=self.fontsize)
    if not nokey:
      #self.leg = axbars.legend(loc='upper right',fontsize=self.fontsize)
      self.leg = axbars.legend(loc='best',fontsize=self.fontsize)
      self.leg.get_frame().set_alpha(0.4)
    #treat the totals differently
    return bot

  def draw_lineplot(self,ax,data,label):
    #ax=axbars.twinx()
    import numpy as np
    ind=np.arange(len(data))+self.barwidth/2.
    ax.plot(ind,data,'o-',color='red',label=label)
    ax.set_ylabel(label,fontsize=self.fontsize)

  def replot(self,label):
    self.vals=label
    #print self.vals
    for i,fig in enumerate(self.figures):
      fx,ax=fig
      cat=ax.get_title()
      ax.cla()
      if i==0:
        data=self.collect_categories(self.scf,label)
      else:
        data=self.find_items(cat,self.scf)
      if data != []: self.draw_barfigure(fx,ax,data,cat)

  def unbalanced(self,val):
    """Criterion for unbalancing"""
    return val > 1.1 or val < 0.9

  def find_unbalanced(self,data):
    """Determine lookup array of unbalanced categories"""
    ipiv=[]
    for i in range(len(data)):
      if self.unbalanced(data[i]):
        ipiv.append(i)
    return ipiv

  def load_unbalancing(self,ax,dict,hosts):
    """Extract the data for plotting the hostname balancings between different categories in bar chart"""
    import numpy as np
    from pylab import cm
    width=0.50
    plts=[]
    key_legend=[]
    values_legend=[]
    icol=1.0
    show_initially=[]
    for cat in self.classes:
      if cat=='Unknown': continue
      try:
        unb=np.array(dict[cat][2:])
        #print 'unbalancing',unb
        #unb2=self.find_unbalanced(unb)
        #print 'unbalanced objects',cat
        #if hosts is not None and (cat=='Convolutions' or cat =='Communications'):
          #print 'unb2',unb2,len(unb),len(hosts)
          #print 'vals',[ [i,unb[i],hosts[i]] for i in unb2]
        ind=np.arange(len(unb))
        pltmp=ax.bar(ind,unb,width,color=cm.jet(icol/len(self.classes)))
        plts.append(pltmp)
        key_legend.append(pltmp[0])
        values_legend.append(cat)
        show_initially.append(dict[cat][0] > 5.0)
        icol+=1.0
        if (width > 0.05):
          width -= 0.05
      except Exception as e:
        print('EXCEPTION FOUND',e)
        print("cat",cat,"not found in workload data")
        return None

    if len(ind) > 2:
      tmp=np.array(hosts) if hosts is not None else None
    else:
      tmp=np.array(["max","min"])
    ax.set_ylabel('Load Unbalancing wrt average')
    ax.set_title('Work Load of different classes')
    if tmp is not None:
      ax.set_xticks(ind+width/2.)
      ax.set_xticklabels(tmp,rotation=90,verticalalignment='bottom')
    ax.set_yticks(np.arange(0,2,0.25))
    ax.axhline(1.0,color='k',linestyle='--')
    leg=ax.legend(np.array(key_legend),np.array(values_legend))#,fancybox=True, shadow=True)
    #leg.get_frame().set_alpha(0.4)
    #show only the most relevant categories
    return leg,plts,show_initially

  def _onpick_curve(self,index,fig,event):
        # on the pick event, find the orig line corresponding to the
    # legend proxy line, and toggle the visibility
    legline = event.artist
    origline = self.lined[index][legline]
    self._toggle_visible(origline,legline)
    fig.canvas.draw()

  def _toggle_visible(self,data,legend):
    for val in data:
      vis = not val.get_visible()
      val.set_visible(vis)
    #last value is like all the others
    if vis:
        legend.set_alpha(1.0)
    else:
        legend.set_alpha(0.2)
