% WHOI letterhead class
%
% Use the letterhead in the following way:
%
% \documentclass[12pt]{WHOIletter}
%
% Then, go to the section below with name, phone number, email,
% etc. and hard-code the new options.

\ProvidesClass{WHOIletter}
\NeedsTeXFormat{LaTeX2e}

\RequirePackage{ifthen}
\RequirePackage{graphics}
\RequirePackage{fancyhdr}

\fancypagestyle{firstpage}{%
  \fancyhf{} % clear all six fields
  \renewcommand{\headrulewidth}{0pt}
  \renewcommand{\footrulewidth}{0pt}
}
\fancypagestyle{followingpage}{%
  \fancyhf{} % clear all six fields
  \renewcommand{\headrulewidth}{0pt}
  \renewcommand{\footrulewidth}{0pt}
}
\pagestyle{followingpage}
\AtBeginDocument{\thispagestyle{firstpage}}

\DeclareOption*{\PassOptionsToClass{\CurrentOption}{letter}}
\ProcessOptions
\LoadClass{letter}

\newcommand*{\WHOIname}[1]{\def\WHOIname{#1}}
\newcommand*{\position}[1]{\def\WHOIposition{#1}}
\newcommand*{\officenum}[1]{\def\WHOIofficenum{#1}}
\newcommand*{\officephone}[1]{\def\WHOIofficephone{#1}}
\newcommand*{\cellphone}[1]{\def\WHOIcellphone{#1}}
\newcommand*{\email}[1]{\def\WHOIemail{#1}}

\def\WHOILogofile{WHOI_newlogo}
\def\WHOIlogo{\includegraphics[width=7.5cm]{\WHOILogofile}}
\def\WHOIheader{
  {\fontfamily{\sfdefault}
    \hspace{-1cm}
    \begin{tabular*}{17.75cm}{@{}r@{}r@{}}
      \parbox{7.8cm}{\resizebox{7.5cm}{!}{\WHOIlogo}}&
      \begin{minipage}[r]{10cm}
        \raggedleft
        \vspace{0.1cm}
        {\large
          \ifx\WHOIname\@empty \else \textbf{\WHOIname}\fi
          \ifthenelse{\isundefined{\WHOIposition}}{}
          {, \textit{\WHOIposition}}
        }\\
        \small{
          \ifx\WHOIemail\@empty \else \WHOIemail \fi\\
          \ifthenelse{\isundefined{\WHOIofficenum}}{}
          {MS \#\WHOIofficenum,}
          266 Woods Hole Road\\
          Woods Hole, MA 02543
        }
      \end{minipage}
    \end{tabular*}
  }
}

\setlength{\topmargin}{-0.5in}
\setlength{\footskip}{0.5in}
\setlength{\textwidth}{7in}
\setlength{\oddsidemargin}{-0.25in}
\setlength{\evensidemargin}{-0.25in}

\renewcommand{\ps@firstpage}
{
   \setlength{\headheight}{0.9in}
   \setlength{\headsep}{0.10in}
   \setlength{\textheight}{8.5in}
   \renewcommand{\@oddhead}{\WHOIheader}
   \renewcommand{\@evenhead}{}
}

\pagestyle{empty}

\renewcommand{\opening}[1]{
  \thispagestyle{firstpage}
  \vspace*{0.4in}
  \par\noindent\toname \\ \toaddress

  \vspace{1ex}

  \noindent #1
}

\newcommand{\sigfig}{%
    \begin{minipage}{\linewidth}
    \raggedright
    \raisebox{0pt}[0pt][0pt]{%
      \includegraphics[width=100px]{Signature.png}}%
    \end{minipage}}
\signature{\sigfig}