OPTION (ENABLE_DOCUMENTATION "Enable build of documentation" OFF)
OPTION (EXCLUDE_DOC_FROM_ALL "Exclude the target doc from all" ON)

configure_file(version.tex.in ${CMAKE_CURRENT_SOURCE_DIR}/version.tex)

IF (ENABLE_DOCUMENTATION)
    project(opal_user_guide NONE)

    include(UseLATEX.cmake)

    MESSAGE(STATUS "Compiling documentation")
    MESSAGE(STATUS "PDFLATEX_COMPILER is '${PDFLATEX_COMPILER}'")
    MESSAGE(STATUS "MAKEINDEX_COMPILER is '${MAKEINDEX_COMPILER}'")
    MESSAGE(STATUS "BIBER_COMPILER is '${BIBER_COMPILER}'")
    MESSAGE(STATUS "")

    IF (EXCLUDE_DOC_FROM_ALL)
        add_latex_document(opal_user_guide.tex
            INPUTS actions.tex autophase.tex benchmarks.tex changelog.tex control.tex conventions.tex dingbat.sty distribution.tex elements.tex errors.tex examples/diagnostic.in examples/matchedDistribution.in examples/matchedDistribution.output examples/multipacting.in examples/OBLA-Gun/OBLA-Gun.in examples/opal-cycl.in examples/opal-cycl-sc.in examples/particlematterinteraction.in examples/tuning.sh fieldemission.tex fieldmaps.tex fieldsolvers.tex figures/Elements/RBend.tex figures/Elements/SBend.tex figures/Multipacting/HB_Fig4.tikz figures/opalt/flowdiagram.tex footer.tex format.tex geometry.tex header.tex html.sty installation.tex introduction.tex lines.tex macros.tex match.tex multipact.tex opalcycl.tex opal-flavours.tex opal-madx.tex opalmap.tex opalt.tex packages.tex partmatter.tex pdfdraftcopy.sty physics.tex quick-installation.tex switches.tex syntax.tex tables.tex track.tex tutorial.tex uml.tex version.tex wakefields.tex
            BIBFILES bibliography.bib
            CONFIGURE macros.tex
            IMAGE_DIRS figures figures/Benchmarks figures/cyclotron figures/Elements figures/Fieldmaps figures/Gun figures/Multipacting  figures/opalt figures/partmatter figures/traveling-wave-structure
            USE_INDEX
            USE_BIBLATEX
            TARGET_NAME doc
            EXCLUDE_FROM_ALL
            )
    ELSE()
        add_latex_document(opal_user_guide.tex
            INPUTS actions.tex autophase.tex benchmarks.tex changelog.tex control.tex conventions.tex dingbat.sty distribution.tex elements.tex errors.tex examples/diagnostic.in examples/matchedDistribution.in examples/matchedDistribution.output examples/multipacting.in examples/OBLA-Gun/OBLA-Gun.in examples/opal-cycl.in examples/opal-cycl-sc.in examples/particlematterinteraction.in examples/tuning.sh fieldemission.tex fieldmaps.tex fieldsolvers.tex figures/Elements/RBend.tex figures/Elements/SBend.tex figures/Multipacting/HB_Fig4.tikz figures/opalt/flowdiagram.tex footer.tex format.tex geometry.tex header.tex html.sty installation.tex introduction.tex lines.tex macros.tex match.tex multipact.tex opalcycl.tex opal-flavours.tex opal-madx.tex opalmap.tex opalt.tex packages.tex partmatter.tex pdfdraftcopy.sty physics.tex quick-installation.tex switches.tex syntax.tex tables.tex track.tex tutorial.tex uml.tex version.tex wakefields.tex
            BIBFILES bibliography.bib
            CONFIGURE macros.tex
            IMAGE_DIRS figures figures/Benchmarks figures/cyclotron figures/Elements figures/Fieldmaps figures/Gun figures/Multipacting  figures/opalt figures/partmatter figures/traveling-wave-structure
            USE_INDEX
            USE_BIBLATEX
            TARGET_NAME doc
            )
    ENDIF()

    add_latex_document(introduction.tex
        INPUTS dingbat.sty footer.tex header.tex html.sty macros.tex packages.tex pdfdraftcopy.sty switches.tex
        BIBFILES bibliography.bib
        CONFIGURE macros.tex
        IMAGE_DIRS figures figures/Gun
        USE_BIBLATEX
        TARGET_NAME doc_introduction
        EXCLUDE_FROM_ALL
        )

    add_latex_document(actions.tex
        INPUTS dingbat.sty footer.tex header.tex html.sty macros.tex packages.tex pdfdraftcopy.sty switches.tex
        BIBFILES bibliography.bib
        CONFIGURE macros.tex
        IMAGE_DIRS figures
        USE_BIBLATEX
        TARGET_NAME doc_action
        EXCLUDE_FROM_ALL
        )

    add_latex_document(autophase.tex
        INPUTS dingbat.sty footer.tex header.tex html.sty macros.tex packages.tex pdfdraftcopy.sty switches.tex
        BIBFILES bibliography.bib
        CONFIGURE macros.tex
        IMAGE_DIRS figures
        USE_BIBLATEX
        TARGET_NAME doc_autophase
        EXCLUDE_FROM_ALL
        )

    add_latex_document(benchmarks.tex
        INPUTS dingbat.sty footer.tex header.tex html.sty macros.tex packages.tex pdfdraftcopy.sty switches.tex
        BIBFILES bibliography.bib
        CONFIGURE macros.tex
        IMAGE_DIRS figures figures/Benchmarks
        USE_BIBLATEX
        TARGET_NAME doc_benchmarks
        EXCLUDE_FROM_ALL
        )

    add_latex_document(changelog.tex
        INPUTS dingbat.sty footer.tex header.tex html.sty macros.tex packages.tex pdfdraftcopy.sty switches.tex
        BIBFILES bibliography.bib
        CONFIGURE macros.tex
        IMAGE_DIRS figures
        USE_BIBLATEX
        TARGET_NAME doc_changelog
        EXCLUDE_FROM_ALL
        )

    add_latex_document(control.tex
        INPUTS dingbat.sty footer.tex header.tex html.sty macros.tex packages.tex pdfdraftcopy.sty switches.tex
        BIBFILES bibliography.bib
        CONFIGURE macros.tex
        IMAGE_DIRS figures
        USE_BIBLATEX
        TARGET_NAME doc_control
        EXCLUDE_FROM_ALL
        )

    add_latex_document(conventions.tex
        INPUTS dingbat.sty footer.tex header.tex html.sty macros.tex packages.tex pdfdraftcopy.sty switches.tex
        BIBFILES bibliography.bib
        CONFIGURE macros.tex
        IMAGE_DIRS figures
        USE_BIBLATEX
        TARGET_NAME doc_conventions
        EXCLUDE_FROM_ALL
        )

    add_latex_document(distribution.tex
        INPUTS dingbat.sty footer.tex header.tex html.sty macros.tex packages.tex pdfdraftcopy.sty switches.tex
        BIBFILES bibliography.bib
        CONFIGURE macros.tex
        IMAGE_DIRS figures
        USE_BIBLATEX
        TARGET_NAME doc_distribution
        EXCLUDE_FROM_ALL
        )

    add_latex_document(elements.tex
        INPUTS dingbat.sty figures/Elements/RBend.tex figures/Elements/SBend.tex footer.tex header.tex html.sty macros.tex packages.tex pdfdraftcopy.sty switches.tex
        BIBFILES bibliography.bib
        CONFIGURE macros.tex
        IMAGE_DIRS figures figures/cyclotron figures/Elements  figures/traveling-wave-structure
        USE_BIBLATEX
        TARGET_NAME doc_elements
        EXCLUDE_FROM_ALL
        )

    add_latex_document(errors.tex
        INPUTS dingbat.sty footer.tex header.tex html.sty macros.tex packages.tex pdfdraftcopy.sty switches.tex
        BIBFILES bibliography.bib
        CONFIGURE macros.tex
        IMAGE_DIRS figures
        USE_BIBLATEX
        TARGET_NAME doc_errors
        EXCLUDE_FROM_ALL
        )

    add_latex_document(fieldemission.tex
        INPUTS dingbat.sty footer.tex header.tex html.sty macros.tex packages.tex pdfdraftcopy.sty switches.tex
        BIBFILES bibliography.bib
        CONFIGURE macros.tex
        IMAGE_DIRS figures
        USE_BIBLATEX
        TARGET_NAME doc_fieldemission
        EXCLUDE_FROM_ALL
        )

    add_latex_document(fieldmaps.tex
        INPUTS dingbat.sty footer.tex header.tex html.sty macros.tex packages.tex pdfdraftcopy.sty switches.tex
        BIBFILES bibliography.bib
        CONFIGURE macros.tex
        IMAGE_DIRS figures figures/Fieldmaps
        USE_BIBLATEX
        TARGET_NAME doc_fieldmaps
        EXCLUDE_FROM_ALL
        )

    add_latex_document(fieldsolvers.tex
        INPUTS dingbat.sty footer.tex header.tex html.sty macros.tex packages.tex pdfdraftcopy.sty switches.tex
        BIBFILES bibliography.bib
        CONFIGURE macros.tex
        IMAGE_DIRS figures
        USE_BIBLATEX
        TARGET_NAME doc_fieldsolver
        EXCLUDE_FROM_ALL
        )

    add_latex_document(format.tex
        INPUTS dingbat.sty footer.tex header.tex html.sty macros.tex packages.tex pdfdraftcopy.sty switches.tex
        BIBFILES bibliography.bib
        CONFIGURE macros.tex
        IMAGE_DIRS figures
        USE_BIBLATEX
        TARGET_NAME doc_format
        EXCLUDE_FROM_ALL
        )

    add_latex_document(geometry.tex
        INPUTS dingbat.sty footer.tex header.tex html.sty macros.tex packages.tex pdfdraftcopy.sty switches.tex
        BIBFILES bibliography.bib
        CONFIGURE macros.tex
        IMAGE_DIRS figures
        USE_BIBLATEX
        TARGET_NAME doc_geometry
        EXCLUDE_FROM_ALL
        )

    add_latex_document(installation.tex
        INPUTS dingbat.sty footer.tex header.tex html.sty macros.tex packages.tex pdfdraftcopy.sty switches.tex
        BIBFILES bibliography.bib
        CONFIGURE macros.tex
        IMAGE_DIRS figures
        USE_BIBLATEX
        TARGET_NAME doc_installation
        EXCLUDE_FROM_ALL
        )

    add_latex_document(lines.tex
        INPUTS dingbat.sty footer.tex header.tex html.sty macros.tex packages.tex pdfdraftcopy.sty switches.tex
        BIBFILES bibliography.bib
        CONFIGURE macros.tex
        IMAGE_DIRS figures
        USE_BIBLATEX
        TARGET_NAME doc_lines
        EXCLUDE_FROM_ALL
        )

    add_latex_document(match.tex
        INPUTS dingbat.sty footer.tex header.tex html.sty macros.tex packages.tex pdfdraftcopy.sty switches.tex
        BIBFILES bibliography.bib
        CONFIGURE macros.tex
        IMAGE_DIRS figures
        USE_BIBLATEX
        TARGET_NAME doc_match
        EXCLUDE_FROM_ALL
        )

    add_latex_document(multipact.tex
        INPUTS dingbat.sty examples/multipacting.in figures/Multipacting/HB_Fig4.tikz footer.tex header.tex html.sty macros.tex packages.tex pdfdraftcopy.sty switches.tex
        BIBFILES bibliography.bib
        CONFIGURE macros.tex
        IMAGE_DIRS figures figures/Multipacting
        USE_BIBLATEX
        TARGET_NAME doc_multipact
        EXCLUDE_FROM_ALL
        )

    add_latex_document(opalcycl.tex
        INPUTS dingbat.sty examples/matchedDistribution.in examples/matchedDistribution.output examples/tuning.sh footer.tex header.tex html.sty macros.tex packages.tex pdfdraftcopy.sty switches.tex
        BIBFILES bibliography.bib
        CONFIGURE macros.tex
        IMAGE_DIRS figures figures/cyclotron
        USE_BIBLATEX
        TARGET_NAME doc_opalcycl
        EXCLUDE_FROM_ALL
        )

    add_latex_document(opalmap.tex
        INPUTS dingbat.sty footer.tex header.tex html.sty macros.tex packages.tex pdfdraftcopy.sty switches.tex
        BIBFILES bibliography.bib
        CONFIGURE macros.tex
        IMAGE_DIRS figures
        USE_BIBLATEX
        TARGET_NAME doc_opal_map
        EXCLUDE_FROM_ALL
        )

    add_latex_document(opal-flavours.tex
        INPUTS dingbat.sty footer.tex header.tex html.sty macros.tex packages.tex pdfdraftcopy.sty switches.tex
        BIBFILES bibliography.bib
        CONFIGURE macros.tex
        IMAGE_DIRS figures figures/cyclotron figures/Fieldmaps
        USE_BIBLATEX
        TARGET_NAME doc_opal_flavours
        EXCLUDE_FROM_ALL
        )

    add_latex_document(opalt.tex
        INPUTS dingbat.sty figures/opalt/flowdiagram.tex footer.tex header.tex html.sty macros.tex packages.tex pdfdraftcopy.sty switches.tex
        BIBFILES bibliography.bib
        CONFIGURE macros.tex
        IMAGE_DIRS figures figures/opalt
        USE_BIBLATEX
        TARGET_NAME doc_opalt
        EXCLUDE_FROM_ALL
        )

    add_latex_document(partmatter.tex
        INPUTS dingbat.sty examples/particlematterinteraction.in footer.tex header.tex html.sty macros.tex packages.tex pdfdraftcopy.sty switches.tex
        BIBFILES bibliography.bib
        CONFIGURE macros.tex
        IMAGE_DIRS figures figures/partmatter
        USE_BIBLATEX
        TARGET_NAME doc_partmatter
        EXCLUDE_FROM_ALL
        )

    add_latex_document(physics.tex
        INPUTS dingbat.sty footer.tex header.tex html.sty macros.tex packages.tex pdfdraftcopy.sty switches.tex
        BIBFILES bibliography.bib
        CONFIGURE macros.tex
        IMAGE_DIRS figures
        USE_BIBLATEX
        TARGET_NAME doc_physics
        EXCLUDE_FROM_ALL
        )

    add_latex_document(quick-installation.tex
        INPUTS dingbat.sty footer.tex header.tex html.sty macros.tex packages.tex pdfdraftcopy.sty switches.tex
        BIBFILES bibliography.bib
        CONFIGURE macros.tex
        IMAGE_DIRS figures
        USE_BIBLATEX
        TARGET_NAME doc_quick_installation
        EXCLUDE_FROM_ALL
        )

    add_latex_document(syntax.tex
        INPUTS dingbat.sty footer.tex header.tex html.sty macros.tex packages.tex pdfdraftcopy.sty switches.tex
        BIBFILES bibliography.bib
        CONFIGURE macros.tex
        IMAGE_DIRS figures
        USE_BIBLATEX
        TARGET_NAME doc_syntax
        EXCLUDE_FROM_ALL
        )

    add_latex_document(tables.tex
        INPUTS dingbat.sty footer.tex header.tex html.sty macros.tex packages.tex pdfdraftcopy.sty switches.tex
        BIBFILES bibliography.bib
        CONFIGURE macros.tex
        IMAGE_DIRS figures
        USE_BIBLATEX
        TARGET_NAME doc_tables
        EXCLUDE_FROM_ALL
        )

    add_latex_document(track.tex
        INPUTS dingbat.sty footer.tex header.tex html.sty macros.tex packages.tex pdfdraftcopy.sty switches.tex
        BIBFILES bibliography.bib
        CONFIGURE macros.tex
        IMAGE_DIRS figures
        USE_BIBLATEX
        TARGET_NAME doc_track
        EXCLUDE_FROM_ALL
        )

    add_latex_document(tutorial.tex
        INPUTS dingbat.sty examples/diagnostic.in examples/OBLA-Gun/OBLA-Gun.in examples/opal-cycl.in examples/opal-cycl-sc.in footer.tex header.tex html.sty macros.tex packages.tex pdfdraftcopy.sty switches.tex
        BIBFILES bibliography.bib
        CONFIGURE macros.tex
        IMAGE_DIRS figures figures/cyclotron figures/Gun
        USE_BIBLATEX
        TARGET_NAME doc_tutorial
        EXCLUDE_FROM_ALL
        )

    add_latex_document(wakefields.tex
        INPUTS dingbat.sty footer.tex header.tex html.sty macros.tex packages.tex pdfdraftcopy.sty switches.tex
        BIBFILES bibliography.bib
        CONFIGURE macros.tex
        IMAGE_DIRS figures
        USE_BIBLATEX
        TARGET_NAME doc_wakefields
        EXCLUDE_FROM_ALL
        )

    add_custom_target(doc_separate DEPENDS doc_introduction doc_autophase doc_benchmarks doc_changelog doc_control doc_conventions doc_distribution doc_elements doc_fieldemission doc_fieldmaps doc_fieldsolver doc_format doc_geometry doc_installation doc_lines doc_multipact doc_opalcycl doc_opalt doc_partmatter doc_physics doc_quick_installation doc_syntax doc_track doc_tutorial doc_wakefields)

    install (FILES opal_user_guide.pdf DESTINATION "${CMAKE_INSTALL_PREFIX}/doc")
ENDIF()
