# Interact with data

<div class='alert alert-info'>
    This section shows how to use the features of the panel dashboard. It is used by the <strong>Single page web app</strong> and the <em>Internal Oceans</em>, <em>Passage Problems</em> and <em>Sub basins</em> steps of the <strong>Multi page web app</strong>.
</div>

<p align="center">
    <img src='img/Panel.png' alt="drawing" width="400">
</p>

Variable to display can be selected at the top of the *left panel* in **Variable** section.

## Colormap

Options on colormap can be set in the *left panel* in **Colormaps** section.

<p align="center">
    <img src='img/colormap_options.png' alt="drawing" width="200">
</p>

Min and max values can be defined. Anything below the minimum value will appear light gray and anything above the maximum value will appear black.

Next, the number of color steps in the colormap can be also set.

<div class='alert alert-info'>
    If set to 1 the colormap is continuous
</div>

Finally a mask can be applied, this colors values below the cutoff to gray and above the cutoff to black. 

<div class='alert alert-info'>
    This is the same as setting the min and max colormap values to the same value
</div>


## Graph tools

<p align="center">
    <img src='img/options.png' alt="drawing" width="300">
</p>

At the top right corner of the *figure panel* are a collections of tools, from left to right these are:

- Box selection
- Lasso selection
- Pan
- Box Zoom
- Mouse Scroll
- Save images
- Refresh / reset all graphs to initial zoom
- Hover (Display information when hovering data with the mouse)

## Selection

As mentioned, in order to change the values of the data, the user can select data in two different ways:

- lasso tool (hand-drawn area)
- box selection (squared area)

<div class='alert alert-info'>
    Long holding either type of selection tool allows the user to change the selection type (union, intersection, ...).
</div>
<p align="center">
    <img src='img/Hold_select.png' alt="drawing" height="200">
    <img src='img/Hold_select_2.png' alt="drawing" height="200">
</p>

<div class='alert alert-info'>
    Selections on multiple graphs are the intersection of the selections.
    This can be useful for when editing specific values in a region using a combination of map + histogram selection.
</div>

<div class='alert alert-info'>
    Prefer Box selection(s) where possible as the calculations for the selections are quicker.
</div>

Then selection can be refined by ticking/unticking boxes in the **Change Values** section in the *left panel*:

<p align="center">
    <img src='img/Select_land_ocean.png' alt="drawing" width="150">
</p>


## Apply new values

<p align="center">
    <img src='img/apply_values.png' alt="drawing" width="200">
</p>

After point(s) have been selected on the graphs the user has the possibility to change the values in three different ways by using buttons in the **Change Values** section in the *left panel*.

If we consider the value in the **Replacement Value** field = `X_value`:

- Absolute (`New_value = X_value`)
- Relative (`New_value = base_value + X_value`)
- Percentage (`New_value = base_value + base_value * X_value %` -> value in **Replacement Value** becomes a percentage)

Apply button can then be clicked to confirm changes. User can also choose to **Undo** or **Redo** changes he has made.

## Save

Once all the changes have been done, the **Save** button can be clicked at the bottom of the *figure panel* to confirm changes and return to the main page.






