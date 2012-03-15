%  PSI-FEL/LEG Note style in LaTeX format                
%
%  Version 1.0  (14/2/95)  JMJ
%  Adapted to latex2e from the old lepnote.sty, just a rather minor
%  modification of 12pt article class.
%
%  Last Modification:13/2/95
%
% Derived from SLS-Note style by Andreas Adelmann
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\def\fileversion{1.0}
\def\filedate{2007/02/13}
\def\docdate {2007/03/13}
\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{sls-note}[\filedate\space sls-note class]
\DeclareOption*{\PassOptionsToClass{\CurrentOption}{article}}
\ProcessOptions

\LoadClass[12pt]{article}
\RequirePackage{a4}
\RequirePackage{graphicx}
\RequirePackage[dvips]{epsfig,rotating}
\RequirePackage{color}
\RequirePackage{latexsym}
\RequirePackage[english]{babel}

\typeout{Document Style `psi-fel/leg-note' 0.1  Febuary 2007}
%it uses the following info
%   \typist{..}  (may be filename) printed at top left of title page
%   \title{..}   
%   \author{..}
%   \date{..}  optional
%   \documentlabel{..}  printed at top right of title page

% You also need to use \begin{document} and \end{document}

 
\def\section{\@startsection {section}{1}{\z@}{-3.5ex plus -1ex minus
    -.2ex}{2.3ex plus .2ex}{\large\bf}}

\def\subsection{\@startsection{subsection}{2}{\z@}{-3.25ex plus -1ex minus
   -.2ex}{1.5ex plus .2ex}{\normalsize\bf}}
\def\subsubsection{\@startsection{subsubsection}{3}{\z@}{-3.25ex plus
 -1ex minus -.2ex}{1.5ex plus .2ex}{\normalsize\bf}}
\def\paragraph{\@startsection
   {paragraph}{4}{\z@}{-2.25ex plus 1ex minus .2ex}{1em}{\normalsize\bf}}
\def\subparagraph{\@startsection
     {subparagraph}{4}{\parindent}{3.25ex plus 1ex minus
     .2ex}{-1em}{\normalsize\bf}}
\setcounter{secnumdepth}{3}
 
% APPENDIX
%
% The \appendix command must do the following:
%    -- reset the section and subsection counters to zero
%    -- redefine the section counter to produce appendix numbers
%    -- redefine the \section command if appendix titles and headings are
%       to look different from section titles and headings.
 
\def\appendix{\par
  \setcounter{section}{0}
  \setcounter{subsection}{0}
  \def\thesection{\Alph{section}}}

%-----------------------------------------------------------------------
% THE REMAINDER OF THIS FILE CONSISTS OF ADDITIONS TO ART12.DOC RATHER
% THAN JUST MODIFICATIONS.
%                           * * Added at 10:07:32 on 3 Apr 1987 by JOWETT
 
 
% Copy following definition from ARTICLE.DOC so that we can set the
% title in \Large rather than \LARGE (which is too large for us).
%
\def\@maketitle{\newpage
%\null
%        * * Above line commented out at 09:37:05 on 4 Dec 1987 by JOWETT
%
% Here we put the typists name at top left and the document  label at the
% top right.  Defaults are an empty box for typist and "Preliminary
% draft" with time of day and date for the label.
% Tabular environments are used to allow for the
% possibility of multiple line document labels - lines separated by \\
% in the argument of \documentlabel - and multiple line typists (which
% could have other uses ...
%
%%                         Modified at 11:19:49 on 26 Sep 1988 by JOWETT
  \begin{tabular}{l} \footnotesize\tt\@typist
  \end{tabular}  \hfill   \begin{tabular}{r} \@documentlabel
                          \end{tabular}
 \vskip 2.5cm               % Vertical space above title.
%                   * * 2em -> 2.5cm at 14:29:08 on 28 Apr 1987 by JOWETT
 \begin{center}
  {\Large\bf \@title \par}     % Title set in \Large, not \LARGE, size.
  \vskip 1.5em                % Vertical space after title.
  {\large                        % each author set in \large, in a
   \lineskip .5em           % tabular environment
   \begin{tabular}[t]{c}\@author
   \end{tabular}\par}
  \vskip 1em              % Vertical space after author.
  {\normalsize \@date}           % Date set in \large size.
\end{center}
 \par
 \vskip 1.5em}                % Vertical space after date.
 
% compute the time in hours and minutes; make new variables \timehh and \timemm
%%                            Added at 11:19:59 on 26 Sep 1988 by JOWETT
%
\newcount\timehh\newcount\timemm
\timehh=\time
\divide\timehh by 60 \timemm=\time
\count255=\timehh\multiply\count255 by -60 \advance\timemm by \count255
%
\def\timeofday{\ifnum\timehh<10 0\fi\number\timehh:%
               \ifnum\timemm<10 0\fi\number\timemm
               }
 
\def\documentlabel#1{\gdef\@documentlabel{#1}}          % see \@maketitle
\gdef\@documentlabel{\tt Preliminary draft\ \timeofday\ \today}
%%                         Modified at 17:59:50 on 23 Sep 1988 by JOWETT
 
\def\typist#1{\gdef\@typist{#1}}          % see \@maketitle
\gdef\@typist{\hbox{}}
% The \documentlabel command does not exist in ordinary LaTeX.
%-----------------------------------------------------------------------

\def\PSI{Paul Scherrer Institut, CH-5232 Villigen PSI, Switzerland}

\def\felabstract#1{\sf Abstract\\[1ex] {\small \sf {#1}}}

\def\feltitle#1#2#3#4#5{

  \setcounter{page}{0}
  \thispagestyle{empty}

  \noindent
  \epsfig{file=PSIlogoBW.pdf, width=95mm}\hspace{0.5cm} \hfill
%\raisebox{11.5mm}{  \epsfig{file=logo-PSIFEL.pdf, width=40mm, angle=-90}}


  \vspace{12mm}
  \begin{flushright}
    {\sf {#1}\\ \today}     % FEL document number
  \end{flushright}

  \begin{center}
    \vfill
    {\LARGE \sf {#2}}\\     % Title
    \vspace{15mm}
    {\Large \sf {#3}}\\     % Author
    \vspace{8mm}
    {\sf {#4}}       % Affiliation
    \vfill
    \parbox{120mm}{\small \sf {#5}}       % text, e.g. an abstract
    \vfill\vfill
  \end{center}

  \newpage
}
 
\endinput
