import matplotlib as mpl
from matplotlib import pyplot as plt
import numpy as np
from astropy.modeling import models, fitting


from dbsp_drp import p200_arm_redux

class ManualTracingGUI:
    """GUI for manually identifying object traces

    Takes in specs_dict:
    specs_dict[target_name] == {
        'spec': spectrum_2d,
        'edges': [all_left_edges, all_right_edges],
        'traces': traces
    }

    How to use:
    Displays one spectrum at a time
    use left and right arrow keys to go back and forth between spectra
    press m to lay down a (m)anual trace under your cursor
    press d to (d)elete a manual trace near your cursor
    """

    def __init__(self, figure, axes, specs_dict):
        self.axes = axes
        self.figure = figure
        self.canvas = figure.canvas

        plt.rcParams['keymap.back'] = ''              # forward / backward keys to enable (Default: left, c, backspace)
        plt.rcParams['keymap.forward'] = ''           # left handed quick navigation (Default: right, v)

        self.canvas.mpl_connect('key_press_event', self.key_press_callback)

        # init needs to load data in!!!
        self.specs_dict = specs_dict
        self.targets = list(specs_dict)#.keys()
        self.spec_index = 0

        self.manual_dict = {targ: [] for targ in self.targets}

        self._plotted_manual_traces = []

        # Initial plot!
        self.plot()
        plt.show()

        self.manual_dict = {k: {'spat_spec': v} for k, v in self.manual_dict.items() if len(v) > 0}

        # Should we try to find FWHM?
        # lets try:
        for target in self.manual_dict:
            targ_dict = self.specs_dict[target]
            spec = targ_dict['spec']
            traces = targ_dict['traces']
            skysub_resid = self.compute_sky_resid(spec)
            
            #g_init = models.Const1D(amplitude=0.)
            
            for trace in self.manual_dict[target]['spat_spec']:
                plt.plot(np.mean(skysub_resid, axis=0), c='gray')
                for trace in traces:
                    plt.axvline(np.median(trace), c='orange')
                
                spat = float(trace.split(":")[0])
                plt.axvline(spat, c='blue')
                
                g_init = models.Gaussian1D(amplitude=1., mean=spat, stddev=1.) + models.Linear1D(slope=0., intercept=0.)
                fit_g = fitting.LevMarLSQFitter()
                xs = np.arange(spec.sciimg.shape[1])
                mask = np.abs(xs - spat) < 20
                g = fit_g(g_init, xs[mask], np.mean(skysub_resid, axis=0)[mask])
                fwhm = g.stddev_0.value * 2 * np.sqrt(2*np.log(2))
                self.manual_dict[target]['fwhm'] = fwhm
                print(f"Fitted FHWM was found to be {fwhm:.2f} pix.")
                
                # maybe modify spat_spec with fitted mean????
                # not great tbh
                # just tell user to be precise
                plt.plot(g(xs), c='cyan')
                plt.xlim(spat - 20, spat + 20)
                plt.ylim(-g.amplitude_0.value, 2*g.amplitude_0.value)
                plt.show()
                user_fwhm = input("If this looks right, hit enter. If the fit was bad, enter your estimate of this object's FWHM in pixels: ")
                if user_fwhm:
                    try:
                        user_fwhm = float(user_fwhm.strip())
                        self.manual_dict[target]['fwhm'] = user_fwhm
                        print(f"Using manually entered FWHM of {user_fwhm:.2f}")
                    except:
                        print(f"Could not cast user-entered FWHM guess {user_fwhm} to a float.")
                        print(f"Keeping fitted fwhm {fwhm}")
        
        for target in self.manual_dict:
            self.manual_dict[target]['needs_std'] = (self.specs_dict[target]['traces'] is None)

    
    @property
    def target(self):
        return self.targets[self.spec_index]
    
    @property
    def spec(self):
        return self.specs_dict[self.target]['spec']
    
    @property
    def edges(self):
        return self.specs_dict[self.target]['edges']
    
    @property
    def traces(self):
        return self.specs_dict[self.target]['traces']
    
    @property
    def mask(self):
        return self.spec.bpmmask == 0
    
    @property
    def sky_resid(self):
        return (self.spec.sciimg - self.spec.skymodel) * np.sqrt(self.spec.ivarmodel) * self.mask
    
    def compute_sky_resid(self, spec):
        return (spec.sciimg - spec.skymodel) * np.sqrt(spec.ivarmodel) * (spec.bpmmask == 0)
    
    def plot(self):
        self.axes.clear()
        self.canvas.set_window_title(self.target)
        
        p200_arm_redux.save_one2dspec(self.axes, self.sky_resid, self.edges, self.traces)
        plt.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=0, hspace=0)

        #self.figure.tight_layout()
        self.canvas.draw()
        
    def plot_manual_traces(self):
        # remove self._plotted_manual_traces
        for plotted_trace in self._plotted_manual_traces:
#            self.axes.lines.remove(plotted_trace)
            plotted_trace.remove()
        self._plotted_manual_traces = []

        # TODO: use pypeit order to draw manual traces
        for trace in self.manual_dict[self.target]:
            spat = float(trace.split(':')[0])
            spec = float(trace.split(':')[1])
            line = self.axes.axvline(spat, c='blue')
            self._plotted_manual_traces.append(line)
            point = self.axes.scatter(spat, spec, c='blue', marker='x')
            self._plotted_manual_traces.append(point)
            #line = self.axes.axvline(spat/440, c='blue')
            #self._plotted_manual_traces.append(line)
        self.canvas.draw()

    def key_press_callback(self, event):
        if not event.inaxes:
            return
        key = event.key

        if key == 'right':
            # display next spectrum
            self.spec_index = (self.spec_index + 1) % len(self.targets)
            self.plot()
            self.plot_manual_traces()
        elif key == 'left':
            # display previous spectrum
            self.spec_index = (self.spec_index - 1) % len(self.targets)
            self.plot()
            self.plot_manual_traces()
        elif key == 'm':
            # save mouse cursor position
            self.manual_dict[self.target].append(f'{event.xdata}:{event.ydata}')
            self.plot_manual_traces()
        elif key == 'd':
            manual_traces = self.manual_dict[self.target]

            if manual_traces:
                best_trace_ix = None
                best_spat = None
                for i in range(len(manual_traces)):
                    this_spat = float(manual_traces[i].split(':')[0])
                    # 5 pixel threshold?
                    if abs(this_spat - event.xdata) < 10:
                        if best_trace_ix is None:
                            best_trace_ix = i
                            best_spat = this_spat
                        else:
                            if (abs(this_spat - event.xdata) < abs(best_spat - event.xdata)):
                                best_trace_ix = i
                                best_spat = this_spat
                
                manual_traces.pop(best_trace_ix)
                self.plot_manual_traces()
