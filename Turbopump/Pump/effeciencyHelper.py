from matplotlib import pyplot as plt, ticker
from PIL import Image
from classDefs import pumpClass
from matplotlib.widgets import TextBox
import os

script_dir = os.path.dirname(os.path.abspath(__file__))
img_path = os.path.join(script_dir, "Ns-Ds Diagram.png")


def effeciency(comp: pumpClass) -> pumpClass:
    """
    Allows the user to define the efficiency by displaying a Balje (Ns-Ds) Diagram with lines plotted on it
    for the input components specific speed and specific diameter. Allows them to apply a multiplication to
    the efficiency as a separate parameter in the pop-up.

    Args:
        comp(pumpClass): A fluid component, either impeller or inducer. Must at least have the following fields:
        Q (m^3/s)
        H(m)
        d_2(m)
        N(rpm)
    """
    ns = comp.N*(comp.Q**.5)/(comp.H**.75)
    ds = comp.D*(comp.H**.25)/(comp.Q**.5)
    img = Image.open(img_path)
    width, height = img.size
    aspect = width / height

    dpi = 100
    # extra space: left margin (0.08) + right margin (0.05) + bottom widget area (1.5in)
    fig_w = width / dpi + 0.13 * (width / dpi)   # add ~13% for left/right margins
    fig_h = height / dpi + 1.5                    # add 1.5in for widgets

    fig = plt.figure(figsize=(fig_w, fig_h))

    # axes rect: [left, bottom, width, height] in figure fraction
    # leave 0.1 on left for y-axis, 0.05 on right for 10000 label, 0.25 at bottom for widgets
    img_h_frac = (height / dpi) / fig_h
    left = 0.10
    right_margin = 0.06
    w = 1.0 - left - right_margin
    
    top_margin = 0.02          # small gap at top
    widget_area = 0.35         # increased from 0.25 — more room for gap + widgets
    img_h_frac_adjusted = 1.0 - top_margin - widget_area
    
    axes_rect = [left, widget_area, w, img_h_frac_adjusted]

    ax1 = fig.add_axes(axes_rect)
    ax1.imshow(img, aspect='auto')
    ax1.set_axis_off()

    # completely independent overlay — no twinning
    ax2 = fig.add_axes(axes_rect, facecolor='none')
    ax2.set_xscale('log')
    ax2.set_yscale('log')
    ax2.set_xlim(1, 10000)
    ax2.set_ylim(0.1, 100)

    ax2.xaxis.set_ticks_position('bottom')
    ax2.xaxis.set_label_position('bottom')
    ax2.yaxis.set_ticks_position('left')
    ax2.yaxis.set_label_position('left')

    ax2.xaxis.set_major_formatter(ticker.ScalarFormatter())
    ax2.yaxis.set_major_formatter(ticker.ScalarFormatter())

    ax2.set_xlabel('Specific Speed (Ns)', labelpad=10)
    ax2.set_ylabel('Specific Diameter (Ds)', labelpad=10)

    # ticks on bottom and left
    ax2.xaxis.set_ticks_position('bottom')
    ax2.xaxis.set_label_position('bottom')
    ax2.yaxis.set_ticks_position('left')
    ax2.yaxis.set_label_position('left')

    ax2.xaxis.set_major_formatter(ticker.ScalarFormatter())
    ax2.yaxis.set_major_formatter(ticker.ScalarFormatter())

    ax2.set_xlabel('Specific Speed (Ns)', labelpad=10)
    ax2.set_ylabel('Specific Diameter (Ds)', labelpad=10)

    ax2.axvline(x=ns, color='violet', linestyle='--', alpha=0.5)
    ax2.axhline(y=ds, color='violet', linestyle='--', alpha=0.5)
    ax2.scatter(ns, ds, color='darkviolet', zorder=5)

    # resize callback — keep image aspect, keep widget area fixed at 1.5in
    def on_resize(event):
        new_w = event.width / fig.dpi
        img_new_h = new_w / aspect
        fig.set_size_inches(new_w, img_new_h + 1.5, forward=False)

    fig.canvas.mpl_connect('resize_event', on_resize)

    # widgets — positioned in the bottom 1.5in regardless of figure size
    ax_box1 = fig.add_axes([0.15, 0.14, 0.25, 0.07])
    ax_box2 = fig.add_axes([0.60, 0.14, 0.25, 0.07])

    ax_msg = fig.add_axes([0.15, 0.04, 0.70, 0.06])

    text_box1 = TextBox(ax_box1, 'Efficiency (Decimal): ', initial='0.5')
    text_box2 = TextBox(ax_box2, 'Efficiency Multiplier: ', initial='1.0')

    ax_msg.set_axis_off()
    error_text = ax_msg.text(0.5, 0.5, '', color='red', ha='center', va='center',
                             transform=ax_msg.transAxes)

    results = {}

    def on_submit(val):
        try:
            results['efficiency'] = float(text_box1.text)
            results['mult'] = float(text_box2.text)
            error_text.set_text('')
            plt.close()
        except ValueError:
            error_text.set_text('Please enter a valid number!')
            fig.canvas.draw()

    text_box1.on_submit(on_submit)
    text_box2.on_submit(on_submit)

    plt.show()

    comp.eta = results['efficiency']*results['mult']
    print("The efficiency is:")
    print(comp.eta* 100, "%")
    return comp
comp = pumpClass()  # or however pumpClass is instantiated
effeciency(comp)
