from matplotlib.axes import Axes
from matplotlib.widgets import AxesWidget

def show_image(imgfile,title=None):
    """
    Show image file using matplotlib imgread. Useful to bypass the
    Jupyter bug for converting a notebook into a pdf file
    """
    import matplotlib.pyplot as plt
    import matplotlib.image as mpimg
    my_image = mpimg.imread(imgfile)
    fig = plt.figure()
    ax1 = fig.add_subplot(1,1,1)
    ax1.set_axis_off()
    ax1.imshow(my_image)
    if title: plt.title(title)
    plt.show()


#interesting slider for vertical matplotlib
#taken from https://stackoverflow.com/questions/25934279/add-a-vertical-slider-with-matplotlib
class VertSlider(AxesWidget):
    """
    A slider representing a floating point range.

    For the slider to remain responsive you must maintain a
    reference to it.

    Attributes:
      *ax*        : the slider :class:`matplotlib.axes.Axes` instance

      *val*       : the current slider value

      *hline*     : a :class:`matplotlib.lines.Line2D` instance
                     representing the initial value of the slider

      *poly*      : A :class:`matplotlib.patches.Polygon` instance
                     which is the slider knob

      *valfmt*    : the format string for formatting the slider text

      *label*     : a :class:`matplotlib.text.Text` instance
                     for the slider label

      *closedmin* : whether the slider is closed on the minimum

      *closedmax* : whether the slider is closed on the maximum

      *slidermin* : another slider - if not *None*, this slider must be
                     greater than *slidermin*

      *slidermax* : another slider - if not *None*, this slider must be
                     less than *slidermax*

      *dragging*  : allow for mouse dragging on slider

    Call :meth:`on_changed` to connect to the slider event
    """
    def __init__(self, ax, label, valmin, valmax, valinit=0.5, valfmt='%1.2f',
                 closedmin=True, closedmax=True, slidermin=None,
                 slidermax=None, dragging=True, **kwargs):
        """
        Create a slider from *valmin* to *valmax* in axes *ax*.

        Additional kwargs are passed on to ``self.poly`` which is the
        :class:`matplotlib.patches.Rectangle` which draws the slider
        knob.  See the :class:`matplotlib.patches.Rectangle` documentation
        valid property names (e.g., *facecolor*, *edgecolor*, *alpha*, ...).

        Parameters
        ----------
        ax : Axes
            The Axes to put the slider in

        label : str
            Slider label

        valmin : float
            The minimum value of the slider

        valmax : float
            The maximum value of the slider

        valinit : float
            The slider initial position

        label : str
            The slider label

        valfmt : str
            Used to format the slider value, fprint format string

        closedmin : bool
            Indicate whether the slider interval is closed on the bottom

        closedmax : bool
            Indicate whether the slider interval is closed on the top

        slidermin : Slider or None
            Do not allow the current slider to have a value less than
            `slidermin`

        slidermax : Slider or None
            Do not allow the current slider to have a value greater than
            `slidermax`


        dragging : bool
            if the slider can be dragged by the mouse

        """
        AxesWidget.__init__(self, ax)

        self.valmin = valmin
        self.valmax = valmax
        self.val = valinit
        self.valinit = valinit
        self.poly = ax.axhspan(valmin, valinit, 0, 1, **kwargs)

        self.hline = ax.axhline(valinit, 0, 1, color='r', lw=1)

        self.valfmt = valfmt
        ax.set_xticks([])
        ax.set_ylim((valmin, valmax))
        ax.set_yticks([])
        ax.set_navigate(False)

        self.connect_event('button_press_event', self._update)
        self.connect_event('button_release_event', self._update)
        if dragging:
            self.connect_event('motion_notify_event', self._update)
        self.label = ax.text(0.5, 1.03, label, transform=ax.transAxes,
                             verticalalignment='center',
                             horizontalalignment='center')

        self.valtext = ax.text(0.5, -0.03, valfmt % valinit,
                               transform=ax.transAxes,
                               verticalalignment='center',
                               horizontalalignment='center')

        self.cnt = 0
        self.observers = {}

        self.closedmin = closedmin
        self.closedmax = closedmax
        self.slidermin = slidermin
        self.slidermax = slidermax
        self.drag_active = False

    def _update(self, event):
        """update the slider position"""
        if self.ignore(event):
            return

        if event.button != 1:
            return

        if event.name == 'button_press_event' and event.inaxes == self.ax:
            self.drag_active = True
            event.canvas.grab_mouse(self.ax)

        if not self.drag_active:
            return

        elif ((event.name == 'button_release_event') or
              (event.name == 'button_press_event' and
               event.inaxes != self.ax)):
            self.drag_active = False
            event.canvas.release_mouse(self.ax)
            return

        val = event.ydata
        if val <= self.valmin:
            if not self.closedmin:
                return
            val = self.valmin
        elif val >= self.valmax:
            if not self.closedmax:
                return
            val = self.valmax

        if self.slidermin is not None and val <= self.slidermin.val:
            if not self.closedmin:
                return
            val = self.slidermin.val

        if self.slidermax is not None and val >= self.slidermax.val:
            if not self.closedmax:
                return
            val = self.slidermax.val

        self.set_val(val)

    def set_val(self, val):
        import six
        xy = self.poly.xy
        xy[1] = 0, val
        xy[2] = 1, val
        self.poly.xy = xy
        self.valtext.set_text(self.valfmt % val)
        if self.drawon:
            self.ax.figure.canvas.draw_idle()
        self.val = val
        if not self.eventson:
            return
        for cid, func in six.iteritems(self.observers):
            func(val)

    def on_changed(self, func):
        """
        When the slider value is changed, call *func* with the new
        slider position

        A connection id is returned which can be used to disconnect
        """
        cid = self.cnt
        self.observers[cid] = func
        self.cnt += 1
        return cid

    def disconnect(self, cid):
        """remove the observer with connection id *cid*"""
        try:
            del self.observers[cid]
        except KeyError:
            pass

    def reset(self):
        """reset the slider to the initial value if needed"""
        if (self.val != self.valinit):
            self.set_val(self.valinit)


def axis_from_data(fig,ax,data):
  "Transform a data tuple into axis coordinates"
  trans=fig.transFigure.inverted()
  ll=trans.transform(ax.transData.transform(data))
  return list(ll)

def data_from_data(fig,dst,src,data):
  "Transform a data tuple of anothe axis in the figure into data of another axis"
  trans=dst.transData.inverted()
  ll=trans.transform(src.transAxes.transform(data))
  return list(ll)

class AxisSet(Axes):
  def create_twin(self):
    self.twin=self.twinx()
  def cla(self):
    if hasattr(self,'twin'): self.twin.cla()
    super(AxisSet,self).cla()
  @classmethod
  def twinify(cls,ax):
    "Include the axis provided as a radix of the set"
    ax.__class__=cls
    ax.create_twin()

class FigureSet():
  """Container for multiple figures.

  Define a container for a plot with the possiblity to switch between simple and gnuplot plotting

  Arguments:
  title: The title of the master figure
  **kwargs: arguments for the axis instance

  """
  def __init__(self,**kwargs):
    import matplotlib.pyplot as plt
    from futile.Utils import kw_pop
    newkw,title=kw_pop('title','',**kwargs)
    self.title=title
    self.figures=[]
    self.showing=False
    self.add(**newkw)

  def __getitem__(self,idx):
    if isinstance(idx,int):
      return self._get_figure(idx)
    else:
      return self._get_figure(figname=idx)

  def _locate(self,figname):
    for i,fig in enumerate(self.figures):
      if figname==fig['Title']:
        return i,fig
    return None

  def exists(self,figname):
    "True if the Figure exists in the Set"
    return self._locate(figname) is not None

  def invoke(self,idx):
    import matplotlib.pyplot as plt
    "Invoke the Figure if it exists. Return None otherwise"
    if isinstance(idx,int):
      if idx >= len(self.figures): return None
      i=idx
    else:
      if not self.exists(idx): return None
      i,fig=self._locate(idx)
    figname=self.figures[i]['TrueTitle']
    return plt.figure(figname)

  def add(self,**kwargs):
    import matplotlib.pyplot as plt
    toadd=str(len(self.figures)+1)
    title=kwargs.get('title','Figure '+toadd if self.title is '' else self.title)
    if self.title != '' and len(self.figures)>0:
      newtitle=title+' ('+self.title+')'
    else:
      newtitle=title
    quitbutton=kwargs.get('QuitButton',False)
    newfig=plt.figure(newtitle)
    axargs=kwargs.get('axes')
    #cfm=plt.get_current_fig_manager()
    totwin=False
    if 'twinaxes' in kwargs: totwin=kwargs.pop('twinaxes')
    if axargs is None:
      newax=newfig.add_subplot(111)
    else:
      rect=axargs.pop('rect')
      newax=plt.axes(rect,**axargs)
    if totwin:
      AxisSet.twinify(newax)
    #connect the home key execpt for the first figure
    from functools import partial
    newfig.canvas.mpl_connect('key_press_event',partial(self._onkey_home,len(self.figures)))
    figdict={'Figure':newfig,'Axes':newax,'Title':title,'TrueTitle':newtitle}
    if quitbutton: figdict.update(self._put_quitbutton())
    self.figures.append(figdict)
    return newfig,newax

  def _onkey_home(self,ifig,event):
    import matplotlib.pyplot as plt
    number=event.key
    #print "pressed",number,'end'
    #return to the main figure if the key "home" is pressed
    if str(number)=='home':
      #print 'raising window'
      self._raisewindow(0)
    elif number == 'q' or number == 'Q':
      if ifig==0:
        self._quitall()
      else:
        figax=self.figures[ifig]
        plt.close(figax['Figure'])

  def _put_quitbutton(self):
    import matplotlib.pyplot as plt
    from matplotlib.widgets import Button
    quitbutton=Button(plt.axes([0.0, 0.0, 0.1, 0.075]), 'Quit')
    quitbutton.on_clicked(self._onclick_quitButton)
    return {'QuitButton':quitbutton}

  def _raisewindow(self,index):
    import matplotlib.pyplot as plt
    if self.invoke(index) is None: return
    cfm=plt.get_current_fig_manager()
    try:
      cfm.window.activateWindow()
      cfm.window.raise_()
    except:
      cfm.window.attributes('-topmost', True)
      #cfm.window.attributes('-topmost', False)
    return cfm

  def _quitall(self):
    import matplotlib.pyplot as plt
    for figax in self.figures:
      plt.close(figax['Figure'])
    self.showing=False
    #self.figures=[] #dereference everything

  def _onclick_quitButton(self,event):
    print ("Good bye!")
    self._quitall()

  def _get_figure(self,index=0,figname=None):
    if figname is not None:
      for fig in self.figures:
        if figname==fig['Title']: break
    else:
      fig=self.figures[index]
    fx=fig["Figure"]
    ax=fig["Axes"]
    return fx,ax

  def show(self,figname=None):
    import matplotlib.pyplot as plt
    if figname is not None and self.exists(figname):
      fx,ax=self._get_figure(figname=figname)
      if not self.showing: self.show()
      fx.show()
    elif figname is None:
      self.showing=True
      plt.show()
