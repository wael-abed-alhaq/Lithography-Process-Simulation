#Lithography Process: Projection Imaging with different masks using FT
#mask 1 (פולס אחד באמצע המסכה )
#mask 2 Vertical pulses (with constant pulse width(w)) and (constant spacing between pulses(space))(פולסים אנכים (ברוחב פולס קבוע) ו-(המרחק בין הפולסים קבוע))
   #The effect of increasing the pulse width (w)(w השפעת הגדלת המפתח )
   #the effect of changing the radios of LPF(aperture of the lens) which leads to pass more orders through LPF   השפעת שינוי הרדיוס של מסנן (מפתח העדשה)  על כמות הסדרים העוברים דרך מפתח העדשה
#mask 3 vertical and horzintal pulses (with constant pulse width and spacing between pulses)פולסים אנכיים ואופקיים (עם רוחב פולס קבוע ומרחק קבוע בין הפולסים)
#mask 4 complex mask with different pulses in both vertical and horizntal spatil מסכה מורכבת עם פולסים שונים בכיוונים אנכיים ואופקיים



import numpy as np
import matplotlib.pylab as plt
import plotly.graph_objects as go

#פרמטרים
NA =0.33 # NA (המפתח הנומרי )




#mask_1
# # Parameters
# N = 1000 # Size of the mask (גודל המסכה)
# low = 498 # Lower boundary for the mask pattern (הגבול התחתון עבור הפולס בתבנית המסכה )
# high = 502 # Upper boundary for the mask pattern (הגבול העליון עבור הפולס בתבנית המסכה)
# space = 20 # Spacing between lines (הרווח בין הפולסים האופקים במסכה)
# limit = 20 # Limit for the inner boundary הגבלה לגבול הפנימי

# # Initialize the 2D mask array (אתחול המטריצה הדו-ממדית ליצירת המסכה)
# d = np.zeros((N, N)) # (N*N)יצירת מטריצת אפסים בגודל

# for i in range(N): # לולאה עבור כל שורה
#     for j in range(N): # לולאה עבור כל עמודה
#         if j >= low and j <= high: #(גבולות הפולס האנכי )
#             if i >= low and i <= high: #(גבולות הפולס האופקי )
#                 d[i, j] = 1 #(הצבת 1 בתאים שעומדים בתנאים לעיל )
# d


# #mask_2
# Parameters
N = 240            # Size of the mask (גודל המסכה -מספר העמודות והשורות במטריצה)

PITCH = 60 # the distance between the starting point of one pulse and the starting point of the next pulse  (המרחק הכולל בין שני פולסים סמוכים בתבנית)
SLIT_PER =50;    # 50% of PITCH
slit_width=round(PITCH*(SLIT_PER)/100);     # slit width (tranparent) [nm] (רוחב של כל פולס )
line_width=round(PITCH*(100-SLIT_PER)/100); # line width (black) [nm](רוחב האזור שבו אין פולסים, הנמדד בין שני פולסים סמוכים)
print('line_width =',line_width)
print('slit_width =',slit_width)
# Initialize the 2D mask array(N*N מטריצת אפסים בגודל )
d = np.zeros((N, N))


# Generate the pulses(יצירת הפולסים)
for i in range(0, N, PITCH):#יצירת מחזרויות של כמה פולסים אנכיים
    d[:, i:i + slit_width] = 1  # Set a vertical pulse of width 'pulse_width')((w)הגדרת פולס אנכי ברוחב מסוים)
d




# #mask_2 (w השפעת הגדלת המפתח )
# # Parameters
# N = 240            # Size of the mask (גודל המסכה -מספר העמודות והשורות במטריצה)
# PITCH = 60 # the distance between the starting point of one pulse and the starting point of the next pulse  (המרחק הכולל בין שני פולסים סמוכים בתבנית)
# SLIT_PER =70;    # 70% of PITCH
# slit_width=round(PITCH*(SLIT_PER)/100);     # slit width (tranparent) [nm] (רוחב של כל פולס )
# line_width=round(PITCH*(100-SLIT_PER)/100); # line width (black) [nm](רוחב האזור שבו אין פולסים, הנמדד בין שני פולסים סמוכים)
# print('line_width =',line_width)
# print('slit_width =',slit_width)
# # Initialize the 2D mask array(N*N מטריצת אפסים בגודל )
# d = np.zeros((N, N))


# # Generate the pulses(יצירת הפולסים)
# for i in range(0, N, PITCH):#יצירת מחזרויות של כמה פולסים אנכיים
#     d[:, i:i + slit_width] = 1  # Set a vertical pulse of width 'pulse_width')((w)הגדרת פולס אנכי ברוחב מסוים)
# d



# #mask_3
# # Parameters
# N = 200            # Size of the mask (גודל המסכה -מספר העמודות והשורות במטריצה)
# PITCH = 50 # the distance between the starting point of one pulse and the starting point of the next pulse  (המרחק הכולל בין שני פולסים סמוכים בתבנית)
# SLIT_PER =50;    # 50% of PITCH
# slit_width=round(PITCH*(SLIT_PER)/100);     # slit width (tranparent) [nm] (רוחב של כל פולס )
# line_width=round(PITCH*(100-SLIT_PER)/100); # line width (black) [nm](רוחב האזור שבו אין פולסים, הנמדד בין שני פולסים סמוכים)
# print('line_width =',line_width)
# print('slit_width =',slit_width)
# # Initialize the 2D mask array(N*N מטריצת אפסים בגודל )
# d = np.zeros((N, N))

# # Generate the pulses
# for i in range(0, N, PITCH):#יצירת מחזרויות של כמה פולסים אנכיים
#     d[:, i:i + slit_width] = 1  # Set a vertical pulse of width 'pulse_width'((w)הגדרת פולס אנכי ברוחב מסוים)

# for j in range (0 , N, PITCH):#יצירת מחזרויות של כמה פולסים אופקים
#    d[j:j+slit_width,:]=1 #Set a horizontal pulse of width 'pulse_width'((w)הגדרת פולס אופקי ברוחב מסוים)

# #כאשר יוצרים מסכה עם פולסים גם בכיוון אופקי וגם בכיוון אנכי, דפוס ההתאבכות במישור פורייה הופך למורכב יותר
# #מה שמוביל לירידה בעוצמה של הסדרים הראשונים

# d



# #mask_4
# # Parameters(פרמטרים)
# N = 200              # Size of the mask (גודל המסכה)
# low = 98             # Lower boundary for the mask pattern (הגבול התחתון עבור הפולס האנכי בתבנית המסכה )
# high = 102           # Upper boundary for the mask pattern (הגבול העליון עבור הפולס האנכי בתבנית המסכה)
# space = 20           # Spacing between lines (הרווח בין הפולסים האופקים במסכה)
# limit = 30           # Limit for the inner boundary  הגבלה לגבול הפנימי של הפולסים האופקים

# # Initialize the 2D mask array (אתחול המטריצה הדו-ממדית ליצירת המסכה)
# d = np.zeros((N, N)) # (N*N)יצירת מטריצת אפסים בגודל

# # Generate the mask(יצירת המסכה)
# for i in range(N):
#     for j in range(N):
#         if low <= j <= high:  # Vertical stripe in the center(יצירת פולס אנכי במרכז המסכה)
#             d[i, j] = 1
#         if i % space == 0 and i > 0:  # Horizontal stripes at intervals(יצירת פולסים אופקים במרווחים קבועים)
#             if limit <= j <= N - limit:
#                 d[i, j] = 1

# # Display the mask(הצגת המסכה)
# d



#הצגת המסכות השונות

# מציג את הפיקסלים ללא החלקה, כך שכל פיקסל מוצג בדיוק בערך המקורי שלו ובגודל חד וברור interpolation='nearest' יצירת גרף חום (הצגת המסכה בצבעים שונים ), הפרמטר
plt.imshow(d, cmap='viridis', interpolation='nearest')

# הוספת בר צבעים  (Colorbar)
plt.colorbar()

#y להפוך את הכיוון של ציר
plt.gca().invert_yaxis()

# הוספת תוויות לצירים
plt.xlabel('x axis')
plt.ylabel('y axis')

#הוספת כותרת לגרף
plt.title('Mask ')
# הצגת הגרף
plt.show()





#הצגת המסכה בתלת מימד

# y-ו x  יצירת וקטורים עבור הצירים
x = np.arange(d.shape[0]) #המכיל ערכים מהאינדקס 0 ועד לסוף העמודות של המטריצה x יצירת וקטור
y = np.arange(d.shape[1]) #המכיל ערכים מהאינדקס 0 ועד לסוף השורות של המטריצה y יצירת וקטור
x, y = np.meshgrid(x, y) #  של הרשת y -ו x יצירת שני מטריצות דו-ממדיות המייצגות את קואורדינטות של
#המטריצות האלה משמשות ליצירת גרפים תלת-ממדיים

#יצירת גרף שטח תלת-ממדי
fig = plt.figure() #חדש שבו נוכל להוסיף תתי-גרפים וגרפים שונים figure יצירת אובייקט
ax = fig.add_subplot(111, projection='3d') # .הוספת תת-גרף תלת-ממדי לדמות
#(מציין תת-גרף יחיד ברשת של 1*1) 111
ax.plot_surface(x, y, d, cmap='viridis') # יצירת הגרף עם צבעים ספציפיים,
ax.view_init(elev=45, azim=45) # קביעת זווית התצוגה של הגרף

# הוספת בר צבעים (Colorbar)
plt.colorbar(ax.plot_surface(x, y, d, cmap='viridis'), ax=ax)  # הוספת בר צבעים לגרף שטח תלת-ממדי


# Set labels (הגדרת תוויות לצירים)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')

#והצגת הגרף הוספת כותרת לגרף
plt.title('3D Surface Plot of d')
plt.show()




#Fourier Transformation in 2-D
d_shifted = np.fft.fftshift(d) #מבצע הזזה של ספקטרום התדרים כדי למרכז את התדרים הנמוכים
dft2d = np.fft.fft2(d_shifted) #מבצע את ההמרה לפורייה על המטריצה שהוזזה.
dft2d_shifted = np.fft.fftshift(dft2d) #מבצע הזזה נוספת כדי למרכז את התדרים הנמוכים בתוצאה של ההמרה לפורייה.
magnitude_spectrum = np.abs(dft2d_shifted) #נרצה לדעת את עוצמת התדרים השונים בתמונה, ולא את ערכי התדרים המורכבים
# magnitude_spectrum = np.abs(np.fft.fftshift(np.fft.fft2(np.fft.fftshift(d))))

# Plot the magnitude spectrum in the frequency domain ( הצגת ספקטרום העוצמה בתחום התדרים)
plt.imshow(np.log1p(magnitude_spectrum), cmap='viridis')
plt.title('Mask diffraction image')
plt.xlabel('spatial frequency in x axis [1/nm]')
#על מנת שיתחיל את הספירה מלמטה ויתקדם למעלה y להתאים את ציר
plt.gca().invert_yaxis()
plt.ylabel('spatial frequency in y axis[1/nm] ')
plt.colorbar(label='Amplitude')
plt.show()



x = np.arange(magnitude_spectrum.shape[0]) #המכיל ערכים מהאינדקס 0 ועד לסוף העמודות של המטריצה x יצירת וקטור
y = np.arange(magnitude_spectrum.shape[1]) #המכיל ערכים מהאינדקס 0 ועד לסוף השורות של המטריצה y יצירת וקטור
x, y = np.meshgrid(x, y)#  של הרשת y -ו x יצירת שני מטריצות דו-ממדיות המייצגות את קואורדינטות של
#המטריצות האלה משמשות ליצירת גרפים תלת-ממדיים


# Create the 3D surface plot using Plotly  # plotly יצירת גרף שטח תלת ממדי באמצעות
fig = go.Figure(data=[go.Surface(z=magnitude_spectrum, x=x, y=y, colorscale= [[0, 'black'],  # Floor color
                                                                              [0.01, 'darkblue'],  # Weak amplitudes
                                                                              [0.1, 'blue'],
                                                                              [0.5, 'yellow'],
                                                                              [1.0, 'red']])])  # Strong amplitudes

# Update layout for better view #עדכון הפריסה של הגרף ע"י הוספת כותרת ו-הגדרת תוויות לצירים ו- שינוי צבבעים ו- התאמת תצוגה
fig.update_layout(
    title='3D Surface Plot of Mask diffraction image',
    scene=dict(
        xaxis_title='spatial frequency in x axis [1/nm]',
        yaxis_title='spatial frequency in y axis [1/nm]',
        zaxis_title='Amplitude',
        xaxis=dict(range=[x.min(), x.max()]),
        yaxis=dict(range=[y.min(), y.max()]),
        zaxis=dict(range=[magnitude_spectrum.min(), magnitude_spectrum.max()])
    )
)

# Show plot # הצגת הגרף
fig.show()



# Define and apply a filter using NA (הגדרת הפילטר)\

#Apperture stop (AS)
filter_radius = round(NA * (N / 2)) # (הגדרת רדיוס הפילטר )
filter_mask = np.ones((N, N))  # Use ones to allow all frequencies (אתחול מסכת הפילטר עם ערכים של 1)
y, x = np.ogrid[:N, :N] # (N*N עבור המטריצה בגודל  y,x יצירת שני וקטורים של האינדקסים לאורך הצירים )
filter_area = (x - N // 2) ** 2 + (y - N // 2) ** 2 <= filter_radius ** 2
filter_mask[~filter_area] = 0  # Set outside the radius to 0 (לחסום את התדרים מחוץ לרדיוס  הפילטר)

# Apply the filter ( השמת המסנן על ספקטרום התדרים של המסכה שלנו ואז הצגת הגרף המתקבל)
# Light passes thru AS (lens)
filtered_spectrum = dft2d_shifted * filter_mask #dft2d_shifted Includes magnitude and phase
# Plot the apperture of the lens
plt.imshow(np.log1p(np.abs(filtered_spectrum)), cmap='viridis')
plt.title('Frequency Spectrum After the LPF')
plt.xlabel('spatial frequency in x axis [1/nm]')
#על מנת שיתחיל את הספירה מלמטה ויתקדם למעלה y להתאים את ציר
plt.gca().invert_yaxis()
plt.ylabel('spatial frequency in y axis[1/nm] ')
plt.colorbar(label='Amplitude')
plt.show()




# Perform the inverse Fourier Transform to see the effect of the filter (ביצוע התמרת פורייה ההפוכה על מנת לראות את ההשפעה של המסנן )
filtered_spectrum_shifted = np.fft.ifftshift(filtered_spectrum)
inverse_filtered_spectrum_filtered = np.fft.ifft2(filtered_spectrum_shifted)
inverse_filtered_spectrum_shifted = np.fft.ifftshift(inverse_filtered_spectrum_filtered)
inverse_image_filtered = np.abs(inverse_filtered_spectrum_shifted)

# יצירת גרף חום תלת-ממדי של התמונה לאחר החלת המסנן והתמרת פורייה הפוכה תוך התאמת מיקום הכותרת ורב הצבעים
fig = go.Figure(data=go.Heatmap(
    z=inverse_image_filtered,
    colorscale='viridis',
    colorbar=dict(title='Filtered Electrical Field Amplitude'),
    hoverongaps=False
))

# Update layout for better view
fig.update_layout(
    title={
        'text': 'Inverse Image After Applying Filter and Inverse Fourier Transform (electrical field of the image)',
        'y': 0.98,  # Vertical position of the title
        'x': 0.48,  # Horizontal position of the title
        'xanchor': 'center',
        'yanchor': 'top'
    },
    xaxis_title='X axis',
    yaxis_title='Y axis',
    yaxis_scaleanchor="x",  # Locks the aspect ratio
    xaxis_constrain="domain",
    yaxis_constrain="domain",
    coloraxis_colorbar=dict(
        x=0.5,  # Shift the colorbar to the right
        y=0.7,  # Center the colorbar vertically
        yanchor="middle",
        lenmode="pixels",
        len=30  # Length of the colorbar
    ),
    margin=dict(t=50)  # Adjust top margin to make room for title (# התאמת השוליים העליונים כדי לפנות מקום לכותרת)
)

fig.update_traces(
     colorbar=dict(
        title='Filtered Electrical Field Amplitude',
        titleside='top',
        x=0.7,  # Shift the colorbar to the right
        y=0.5,  # Center the colorbar vertically
        len=0.75,  # Adjust the length of the colorbar
        thickness=15,  # Adjust the thickness of the colorbar
    )
)

# Show the plot
fig.show()





# Square the magnitude of the electrical field to get the image intensity-energy(עוצמת התמונה- האנרגיה)
intensity_image = np.square(inverse_image_filtered)


# Create a heatmap with Plotly
fig = go.Figure(data=go.Heatmap(
    z=intensity_image,
    colorscale='viridis',
    colorbar=dict(title='Intensity'),
    hoverongaps=False
))

# Update layout for better view
fig.update_layout(
     title={
        'text': '2D Intensity Image After Applying Filter and Inverse Fourier Transform',
        'y': 0.98,
        'x': 0.48,
        'xanchor': 'center',
        'yanchor': 'top'
    },
    xaxis_title='X axis',
    yaxis_title='Y axis',
    yaxis_scaleanchor="x",  # Locks the aspect ratio(keeping the ratio of width to height consistent in a plot or image.) # שמירה על עקביות היחס בין רוחב לגובה בעלילה או תמונה
    xaxis_constrain="domain",
    yaxis_constrain="domain",
    margin=dict(t=50)  # Adjust top margin to make room for title (# התאמת השוליים העליונים כדי לפנות מקום לכותרת)
)

#adjusting the appearance and position of the colorbar.
fig.update_traces(
     colorbar=dict(
        title='Intensity',
        titleside='top',
        x=0.7,  # Shift the colorbar to the right
        y=0.5,  # Center the colorbar vertically
        len=0.75,  # Adjust the length of the colorbar
        thickness=15,  # Adjust the thickness of the colorbar
    )
)
# Show the plot
fig.show()





# 3D Surface Plot of the Image
x = np.arange(inverse_image_filtered.shape[1])
y = np.arange(inverse_image_filtered.shape[0])
x, y = np.meshgrid(x, y)


# Create the 3D surface plot using Plotly
fig = go.Figure(data=[go.Surface(z=inverse_image_filtered, x=x, y=y, colorscale= [[0, 'black'],  # Floor color
                                                                              [0.01, 'darkblue'],  # Weak amplitudes
                                                                              [0.1, 'blue'],
                                                                              [0.5, 'yellow'],
                                                                              [1.0, 'red']])])  # Strong amplitudes


# Update layout for better view
fig.update_layout(
     title={
        'text': '3D Surface Plot of Intensity Image',
        'y': 0.98,
        'x': 0.48,
        'xanchor': 'center',
        'yanchor': 'top'
    },
    xaxis_title='X axis',
    yaxis_title='Y axis',
    yaxis_scaleanchor="x",  # Locks the aspect ratio(keeping the ratio of width to height consistent in a plot or image.) # שמירה על עקביות היחס בין רוחב לגובה בעלילה או תמונה
    xaxis_constrain="domain",
    yaxis_constrain="domain",
    margin=dict(t=50)  # Adjust top margin to make room for title (# התאמת השוליים העליונים כדי לפנות מקום לכותרת)
)



fig.update_traces(
     colorbar=dict(
        title='Intensity',
        titleside='top',
        x=0.7,  # Shift the colorbar to the right
        y=0.5,  # Center the colorbar vertically
        len=0.75,  # Adjust the length of the colorbar
        thickness=15,  # Adjust the thickness of the colorbar
    )
)
# Show the plot
fig.show()

